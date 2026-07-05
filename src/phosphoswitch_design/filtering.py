"""
filtering.py — plausibility filter, H1/H2 gate, and diversity selection.

Stage 02 plausibility filter
-----------------------------
Rejects sequences that are biologically implausible before any expensive
computation:

    - Length must be exactly 59 aa
    - Position 30 must be Y (the phospho-Tyr; kinase recognition)
    - Position 59 must be V (C-terminal anchor residue)
    - At most MAX_MUTATIONS vs WT template
    - Composition: hydrophobic %, leucine %, charged %, aromatic %
    - No homopolymer run of 5+ identical residues
    - Capped K+R count, D+E count, total charged count

These are intentionally set generous relative to WT (LMNA construct has 22
total charged residues; default MAX_TOTAL_CHARGED is 26) so the filter does
not falsely reject interesting designs, just obvious garbage.

WT composition reference (59 aa): ASSTPLSPTRITRLQEKEDLQELNDRLAVYIDRVRSLETENAGLRLRITESEEVVSREV
  Hydrophobic (LIVMFAW): 16 → 27%
  Leucine:                7 → 12%
  Charged (DEKR):        22 → 37%
  Aromatic (FWY):         2 →  3%

Stage 03 diversity selection
-----------------------------
After the plausibility + mechanism scoring in stage 02, stage 03 selects a
manageable subset (~10,000) for deep Rosetta validation, using:

  1. H1/H2 directional contact filter — eliminate designs whose mechanism
     contradicts design intent (this is the most important gate)
  2. Rank by |deviation_score_diff| within each direction pool
  3. Diversity caps: max N per subspace, max N per mutation signature

Direction allocation: 40% hairpin-favoring, 40% straight-favoring, 20% neutral.
"""

from __future__ import annotations
from collections import defaultdict

# ---------------------------------------------------------------------------
# LMNA construct constants
# ---------------------------------------------------------------------------
TEMPLATE = "ASSTPLSPTRITRLQEKEDLQELNDRLAVYIDRVRSLETENAGLRLRITESEEVVSREV"
EXPECTED_LEN = 59

# Composition sets
HYDROPHOBIC = frozenset("LIVMFAW")
CHARGED     = frozenset("DEKR")
AROMATIC    = frozenset("FWY")

# Default plausibility thresholds
# WT reference: 22 charged, 10 K+R, 12 D+E, 27% hydrophobic, 12% leu, 3% aro
MAX_MUTATIONS      = 12
MAX_HYDROPHOBIC_PCT = 0.50
MAX_LEUCINE_PCT    = 0.18
MIN_CHARGED_PCT    = 0.18
MAX_AROMATIC_PCT   = 0.12
MAX_RUN_OF_SAME    = 4
MAX_TOTAL_CHARGED  = 26   # WT=22, generous ceiling
MAX_TOTAL_KR       = 14   # WT=10
MAX_TOTAL_DE       = 14   # WT=12


# ---------------------------------------------------------------------------
# Plausibility filter
# ---------------------------------------------------------------------------
def has_homopolymer_run(seq: str, n: int = MAX_RUN_OF_SAME) -> bool:
    """Return True if *seq* contains a run of N+ identical residues."""
    for i in range(len(seq) - n + 1):
        if len(set(seq[i:i + n])) == 1:
            return True
    return False


def passes_plausibility(
    seq: str,
    max_mutations: int = MAX_MUTATIONS,
    max_hydrophobic_pct: float = MAX_HYDROPHOBIC_PCT,
    max_leucine_pct: float = MAX_LEUCINE_PCT,
    min_charged_pct: float = MIN_CHARGED_PCT,
    max_aromatic_pct: float = MAX_AROMATIC_PCT,
    max_run_of_same: int = MAX_RUN_OF_SAME,
    max_total_charged: int = MAX_TOTAL_CHARGED,
    max_total_kr: int = MAX_TOTAL_KR,
    max_total_de: int = MAX_TOTAL_DE,
    relaxed: bool = False,
) -> tuple[bool, str]:
    """Apply all plausibility checks.

    Parameters
    ----------
    seq : str
        Candidate amino acid sequence.
    relaxed : bool
        If True, skip composition percentage and charged-count filters.
        Still enforces: length=59, Y30, V59, max_mutations, no homopolymer≥5.
        Use for exploratory analyses where wet-lab assessment replaces these gates.

    Returns
    -------
    (passed, reason) where *reason* is 'ok' or a short code for rejection.
    """
    if len(seq) != EXPECTED_LEN:
        return False, f"len={len(seq)}"
    if seq[29] != 'Y':
        return False, "pos30_not_Y"
    if seq[58] != 'V':
        return False, "pos59_not_V"

    n_muts = sum(1 for t, s in zip(TEMPLATE, seq) if t != s)
    if n_muts > max_mutations:
        return False, f"muts={n_muts}"

    if has_homopolymer_run(seq, max_run_of_same + 1):
        return False, "homopolymer"

    if not relaxed:
        L = len(seq)
        n_hyd = sum(1 for c in seq if c in HYDROPHOBIC)
        n_leu = seq.count('L')
        n_chg = sum(1 for c in seq if c in CHARGED)
        n_aro = sum(1 for c in seq if c in AROMATIC)

        if n_hyd / L > max_hydrophobic_pct:
            return False, f"hyd={n_hyd/L:.2f}"
        if n_leu / L > max_leucine_pct:
            return False, f"leu={n_leu/L:.2f}"
        if n_chg / L < min_charged_pct:
            return False, f"chg={n_chg/L:.2f}"
        if n_aro / L > max_aromatic_pct:
            return False, f"aro={n_aro/L:.2f}"

        n_KR = seq.count('K') + seq.count('R')
        n_DE = seq.count('D') + seq.count('E')
        n_charged_total = n_KR + n_DE

        if n_charged_total > max_total_charged:
            return False, f"chg={n_charged_total}"
        if n_KR > max_total_kr:
            return False, f"KR={n_KR}"
        if n_DE > max_total_de:
            return False, f"DE={n_DE}"

    return True, "ok"


# ---------------------------------------------------------------------------
# Mutation signature (for diversity grouping)
# ---------------------------------------------------------------------------
def get_mutation_signature(seq: str, template: str = TEMPLATE) -> str:
    """Return a region-based mutation pattern string for diversity grouping.

    Four regions:
        nterm     (1-14)   N-terminal helix
        phos_loop (15-30)  Phospho-Tyr loop
        central   (31-43)  Central helix
        tail      (44-59)  C-terminal hairpin face

    Example: 'n0_p2_c1_t3' means 0 N-term mutations, 2 phos-loop mutations,
    1 central mutation, 3 tail mutations.
    """
    if seq == template:
        return "WT"
    muts = [
        (i + 1, t, s)
        for i, (t, s) in enumerate(zip(template, seq))
        if t != s
    ]
    if not muts:
        return "WT"

    regions = {'nterm': 0, 'phos_loop': 0, 'central': 0, 'tail': 0}
    for p, _, _ in muts:
        if 1 <= p <= 14:
            regions['nterm'] += 1
        elif 15 <= p <= 30:
            regions['phos_loop'] += 1
        elif 31 <= p <= 43:
            regions['central'] += 1
        elif 44 <= p <= 59:
            regions['tail'] += 1

    return (f"n{regions['nterm']}_p{regions['phos_loop']}"
            f"_c{regions['central']}_t{regions['tail']}")


# ---------------------------------------------------------------------------
# H1/H2 directional contact filter
# ---------------------------------------------------------------------------
def passes_h1_h2_filter(record: dict) -> bool:
    """Return True if the candidate's phosphate contacts match its hypothesis.

    H1 (phospho stabilises STRAIGHT):
        HP donor count decreased from WT  OR  ST donor count increased from WT.

    H2 (phospho stabilises HAIRPIN):
        HP donor count increased from WT  OR  ST donor count decreased from WT.

    Parameters
    ----------
    record : dict
        A row dict containing at least:
        - 'hypothesis': 'H1' or 'H2'
        - 'hp_donors_changed': int  (HP donors this design − HP donors WT)
        - 'st_donors_changed': int  (ST donors this design − ST donors WT)
    """
    hyp = record.get('hypothesis', '?')
    hp_change = float(record.get('hp_donors_changed', 0))
    st_change = float(record.get('st_donors_changed', 0))

    if hyp == 'H1':
        return (hp_change < 0) or (st_change > 0)
    if hyp == 'H2':
        return (hp_change > 0) or (st_change < 0)
    return False


# ---------------------------------------------------------------------------
# Diversity-capped selection
# ---------------------------------------------------------------------------
def select_with_diversity(
    pool: list[dict],
    n_target: int,
    label: str,
    max_per_signature: int = 30,
    max_per_subspace: int = 1000,
) -> tuple[list[dict], dict, dict]:
    """Greedily select up to *n_target* records from *pool* with diversity caps.

    The pool must already be sorted in preference order (highest deviation
    first).  A record is accepted if and only if both its mutation signature
    count and its subspace count are below their caps.

    Parameters
    ----------
    pool : list[dict]
        Pre-sorted candidate records.
    n_target : int
        Maximum number to select.
    label : str
        Tag written into each selected record's '_selection_pool' field.
    max_per_signature : int
        Maximum accepted records sharing the same mutation signature.
    max_per_subspace : int
        Maximum accepted records from the same subspace.

    Returns
    -------
    (selected, sig_counts, sub_counts)
    """
    selected: list[dict] = []
    sig_counts: dict[str, int] = defaultdict(int)
    sub_counts: dict[str, int] = defaultdict(int)

    for r in pool:
        sig = get_mutation_signature(r.get('sequence', ''))
        sub = r.get('subspace', '?')

        if sig_counts[sig] >= max_per_signature:
            continue
        if sub_counts[sub] >= max_per_subspace:
            continue

        r['_selection_pool'] = label
        selected.append(r)
        sig_counts[sig] += 1
        sub_counts[sub] += 1

        if len(selected) >= n_target:
            break

    return selected, dict(sig_counts), dict(sub_counts)
