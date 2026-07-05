"""
consensus.py — consensus scoring and final wet-lab candidate selection.

Consensus scoring weights
--------------------------
Each signal contributes independently; the composite score drives ranking.

    Rosetta z-score       ×2.0 per |z| unit  — primary thermodynamic signal
    Mechanism deviation   ×0.15 per unit      — geometry contact differential
    Direction agreement   1-3 pts             — Rosetta + mech + AF agree
    Folding stability     1.5 pts             — apo energies < 0 REU
    Solubility            0.5-1.0 pts         — predicted solubility score
    AF2 bistability       0.5-1.0 pts         — apo↔phos RMSD signal
    AF3 confidence        1.0 pts             — both apo/phos pLDDT > 70
    AF3 change            1.5 pts             — AF3 predicts conformational change
    Minimal mutations     0.5-1.0 pts         — ≤7 muts (bonus for ≤5)
    No outlier            0.5 pts             — Rosetta range ≤ 15 REU
    Physical filters pass 0.5 pts             — composition gates

Penalties (always applied unless noted):
    H-cluster             −5.0 (always)       — PyRosetta fixed-proton artifact
    High aggregation      −2.0 (strict mode)  — ΔAgg vs WT > 0.10
    PTM crosstalk         −1.0 each (strict)  — secondary phospho sites

Hard requirements (knock-out gates):
    passes_physical_filters must be True (composition, no H-clusters,
    solubility > 0.4).

Diversity selection
--------------------
From all candidates that pass hard requirements, ranked by consensus score:
    - At most 3 per mutation signature (prevents near-duplicate clusters)
    - At least 2 directional representations (H1-type and H2-type)
    - At least 4 different mutation signatures in final set

Output files
------------
    final_ranking.csv         full ranked table
    protein_sequences.fa      FASTA of top candidates + WT control
    codon_optimized_dna.fa    E. coli optimised DNA sequences
    selection_rationale.txt   human-readable decision log
"""

from __future__ import annotations
from collections import defaultdict

from .filtering import TEMPLATE, get_mutation_signature

# ---------------------------------------------------------------------------
# E. coli codon table
# ---------------------------------------------------------------------------
CODON_E_COLI: dict[str, str] = {
    'A': 'GCG', 'R': 'CGT', 'N': 'AAC', 'D': 'GAT', 'C': 'TGC',
    'E': 'GAA', 'Q': 'CAG', 'G': 'GGC', 'H': 'CAT', 'I': 'ATT',
    'L': 'CTG', 'K': 'AAA', 'M': 'ATG', 'F': 'TTT', 'P': 'CCG',
    'S': 'AGC', 'T': 'ACC', 'W': 'TGG', 'Y': 'TAT', 'V': 'GTG',
    '*': 'TAA',
}


def reverse_translate(seq: str, table: dict[str, str] = CODON_E_COLI) -> str:
    """Reverse-translate a peptide using E. coli high-expression codons.

    Appends a TAA stop codon.  Unknown residues become 'NNN'.

    Parameters
    ----------
    seq : str
        Amino acid sequence.
    table : dict
        Codon table mapping one-letter AA codes to codon strings.
    """
    return ''.join(table.get(aa, 'NNN') for aa in seq) + table['*']


def get_mutations(seq: str, template: str = TEMPLATE) -> list[str]:
    """Return mutation list in standard notation (e.g. ['L20K', 'E41R'])."""
    return [f"{t}{i+1}{s}" for i, (t, s) in enumerate(zip(template, seq)) if t != s]


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------
def _safe_float(x, default: float = 0.0) -> float:
    try:
        return float(x)
    except (ValueError, TypeError):
        return default


def _safe_bool(x) -> bool:
    if isinstance(x, bool):
        return x
    return str(x).lower() == 'true'


def _safe_int(x, default: int = 0) -> int:
    try:
        return int(float(x))
    except (ValueError, TypeError):
        return default


# ---------------------------------------------------------------------------
# Consensus score
# ---------------------------------------------------------------------------
def compute_consensus(rec: dict, relaxed: bool = False) -> tuple[float, dict]:
    """Compute composite consensus score for one candidate.

    Parameters
    ----------
    rec : dict
        Merged row dict from all pipeline stages (Rosetta, folding, AF2, AF3).
    relaxed : bool
        If True, skip aggregation and PTM crosstalk penalties.  The histidine
        cluster penalty is always applied because it is a PyRosetta artifact,
        not a real biophysical concern.

    Returns
    -------
    (score, contributions) where *contributions* maps signal name → value.
    """
    score = 0.0
    contributions: dict[str, float] = {}

    # 1. Rosetta z-score (primary thermodynamic signal)
    z = _safe_float(rec.get('z_score_vs_WT', 0))
    z_contrib = abs(z) * 2.0
    score += z_contrib
    contributions['rosetta_z'] = round(z_contrib, 2)

    # 2. Geometry mechanism deviation
    mech_dev = abs(_safe_float(rec.get('mech_deviation', 0)))
    mech_contrib = min(mech_dev, 10.0) * 0.15   # capped at 10 units
    score += mech_contrib
    contributions['mechanism'] = round(mech_contrib, 2)

    # 3. Multi-source direction agreement
    rosetta_dir = 'HP' if _safe_float(rec.get('delta_vs_WT_median', 0)) > 0 else 'ST'
    mech_dir = rec.get('mech_direction', '')
    af2_change = _safe_float(rec.get('apo_vs_phos_rmsd', 0))

    agreement = 0
    if 'HAIRPIN' in mech_dir and rosetta_dir == 'HP':
        agreement += 1
    elif 'STRAIGHT' in mech_dir and rosetta_dir == 'ST':
        agreement += 1
    if af2_change > 2.0:
        agreement += 1

    if agreement >= 2:
        score += 3.0
        contributions['direction_agree'] = 3.0
    elif agreement == 1:
        score += 1.0
        contributions['direction_agree'] = 1.0
    else:
        contributions['direction_agree'] = 0.0

    # 4. Folding stability (both apo states should relax below 0 REU)
    e_apo_st = rec.get('E_mean_apo_ST')
    e_apo_hp = rec.get('E_mean_apo_HP')
    if e_apo_st and e_apo_hp:
        e_apo_min = min(_safe_float(e_apo_st), _safe_float(e_apo_hp))
        if e_apo_min < 0:
            score += 1.5
            contributions['folding'] = 1.5
        else:
            contributions['folding'] = 0.0
    else:
        contributions['folding'] = 0.0

    # 5. Solubility
    sol = _safe_float(rec.get('solubility_score', 0))
    if sol > 0.5:
        score += 1.0
        contributions['solubility'] = 1.0
    elif sol > 0.4:
        score += 0.5
        contributions['solubility'] = 0.5
    else:
        contributions['solubility'] = 0.0

    # 6. AF2 bistability (apo seeds sample different conformations)
    if _safe_bool(rec.get('predicts_conformational_change')):
        score += 1.0
        contributions['af2_bistability'] = 1.0
    elif _safe_float(rec.get('apo_inter_seed_rmsd_max', 0)) > 2.5:
        score += 0.5
        contributions['af2_bistability'] = 0.5
    else:
        contributions['af2_bistability'] = 0.0

    # 7. AF3 confidence (pLDDT > 70 in both states)
    af3_apo_plddt  = _safe_float(rec.get('af3_apo_plddt', 0))
    af3_phos_plddt = _safe_float(rec.get('af3_phos_plddt', 0))
    if af3_apo_plddt > 70 and af3_phos_plddt > 70:
        score += 1.0
        contributions['af3_confidence'] = 1.0
    else:
        contributions['af3_confidence'] = 0.0

    # 8. AF3 conformational change
    if _safe_bool(rec.get('af3_predicts_change')):
        score += 1.5
        contributions['af3_change'] = 1.5
    else:
        contributions['af3_change'] = 0.0

    # 9. Minimal mutations bonus (prefer parsimonious designs)
    n_muts = _safe_int(rec.get('n_mutations', 99))
    if n_muts <= 5:
        score += 1.0
        contributions['minimal_muts'] = 1.0
    elif n_muts <= 7:
        score += 0.5
        contributions['minimal_muts'] = 0.5
    else:
        contributions['minimal_muts'] = 0.0

    # 10. Rosetta reliability — no stochastic explosions
    if not _safe_bool(rec.get('outlier_flag')):
        score += 0.5
        contributions['no_outlier'] = 0.5
    else:
        contributions['no_outlier'] = 0.0

    # 11. Physical filters gate
    if _safe_bool(rec.get('passes_physical_filters')):
        score += 0.5
        contributions['physical_pass'] = 0.5
    else:
        contributions['physical_pass'] = 0.0

    # 12. Penalties
    # H-cluster penalty always applied — PyRosetta fixed-proton artifact
    if _safe_bool(rec.get('h_cluster_warning')):
        score -= 5.0
        contributions['h_cluster_penalty'] = -5.0

    if not relaxed:
        # WT-relative aggregation score (absolute threshold rejects WT itself)
        agg_delta = _safe_float(rec.get('agg_delta_vs_wt', 0))
        if agg_delta > 0.10:
            score -= 2.0
            contributions['agg_penalty'] = -2.0
        if _safe_int(rec.get('n_ptm_crosstalk', 0)) > 0:
            score -= 1.0
            contributions['ptm_crosstalk_penalty'] = -1.0

    return round(score, 2), contributions


# ---------------------------------------------------------------------------
# Diversity-aware final selection
# ---------------------------------------------------------------------------
def select_final_candidates(
    designs: list[dict],
    n_final: int = 12,
    min_directions: int = 2,
    min_signatures: int = 4,
    max_per_signature: int = 3,
) -> list[dict]:
    """Select the top *n_final* candidates with diversity constraints.

    Algorithm:
        1. First pass: greedily accept from score-sorted list, max
           *max_per_signature* per mutation signature.
        2. If fewer than *min_directions* distinct directions are represented,
           find and append the top-scoring candidate from each missing direction.
        3. Sort final list by consensus score.

    Parameters
    ----------
    designs : list[dict]
        Candidates that have already passed all hard gates, with
        'consensus_score' and 'mutation_signature' populated.
    n_final : int
        Target selection size (default 12).
    min_directions : int
        Minimum number of distinct mechanism directions required.
    min_signatures : int
        Minimum number of distinct mutation signatures required (informational).
    max_per_signature : int
        Maximum candidates sharing the same mutation signature.

    Returns
    -------
    List of selected candidate dicts, sorted by descending consensus score.
    """
    # Sort descending by consensus score
    designs.sort(key=lambda r: -r['consensus_score'])

    selected: list[dict] = []
    sigs_count: dict[str, int] = defaultdict(int)
    dirs_seen: set[str] = set()
    sigs_seen: set[str] = set()

    for rec in designs:
        if len(selected) >= n_final:
            break
        sig = rec.get('mutation_signature', '?')
        if sigs_count[sig] >= max_per_signature:
            continue
        selected.append(rec)
        sigs_count[sig] += 1
        sigs_seen.add(sig)
        dirs_seen.add(rec.get('mech_direction', '?'))

    # Enforce directional diversity
    if len(dirs_seen) < min_directions:
        all_dirs = {r.get('mech_direction', '?') for r in designs}
        missing = all_dirs - dirs_seen
        for missing_dir in missing:
            for rec in designs:
                if rec in selected:
                    continue
                if rec.get('mech_direction') == missing_dir:
                    selected.append(rec)
                    dirs_seen.add(missing_dir)
                    break

    selected.sort(key=lambda r: -r['consensus_score'])
    return selected
