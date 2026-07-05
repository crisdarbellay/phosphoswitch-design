"""
mechanism.py — phosphate contact scoring and H1/H2 hypothesis classification.

Biology context
---------------
LMNA Y45 is phosphorylated by Src kinase.  The ±29-aa construct (59 aa total)
has Y45 at position 30.  State A is the straight extended-helix backbone;
State B is the pulled hairpin backbone.

Two competing hypotheses about what phosphorylation does:

    H1: phospho stabilises STRAIGHT (helix)
        → phospho must make better contacts on straight backbone than hairpin
        → design target: break/reduce hairpin contacts, enhance straight contacts

    H2: phospho stabilises HAIRPIN
        → phospho must make better contacts on hairpin backbone than straight
        → design target: create hairpin contacts, reduce straight contacts

Scoring approach
----------------
Purely geometry-based — no force-field, no docking.  The phosphate centroid
is located from real PDB coordinates (P/OP atoms of the phospho-Tyr at pos 30).
For every other residue in the sequence, the CB (or CA for Gly) distance to
that centroid is looked up; a per-AA cutoff and weight table converts distances
to a score.

This approach is anti-gameable: it cannot be fooled by stacking copies of the
same AA because it uses real backbone geometry from pre-computed PDBs.

Functions
---------
parse_pdb           Extract phosphate centroid and CA/CB coords from a PDB.
dist                Euclidean distance.
aa_score            Per-residue score at a given distance.
count_interactions  Full site analysis: donors, h-bonds, repellers, regional.
bidirectional_score Score switch magnitude in BOTH directions simultaneously.
mechanism_score     CSV-ready dict with all HP/ST breakdown fields.
classify_h1_h2      Boolean: does this record match its design hypothesis?
"""

from __future__ import annotations
from typing import Optional

# ---------------------------------------------------------------------------
# Phosphate interaction parameters
# Per-AA distance cutoff and binding weight, calibrated to sidechain reach:
#   K  cutoff 10.0 Å because NZ extends ~6.4 Å from CB
#   R  cutoff 11.0 Å because CZ/NH1/NH2 extend ~7.5 Å
#   H  cutoff  9.0 Å because ND1/NE2 extend ~5.7 Å
#   S/T cutoff 6.0 Å because OG/OG1 extend ~2.5 Å
#   N/Q cutoff 8/9 Å  — amide oxygens / nitrogens
#   Y  cutoff 10.0 Å  — phenolic OH at tip of ring
# ---------------------------------------------------------------------------
PHOS_BINDERS: dict[str, tuple[float, float]] = {
    'K': (10.0, 1.5),   # lysine  — electrostatic donor
    'R': (11.0, 2.0),   # arginine — bidentate donor (highest weight)
    'H': (9.0,  0.8),   # histidine — neutral at pH 7, partial donor
    'S': (6.0,  0.5),   # serine  — H-bond donor (short range only)
    'T': (6.0,  0.5),   # threonine
    'Y': (10.0, 0.5),   # tyrosine phenolic OH
    'N': (8.0,  0.3),   # asparagine amide
    'Q': (9.0,  0.3),   # glutamine amide
}

# Negatively charged residues repel the phosphate group (charge-charge clash).
PHOS_REPELLERS: dict[str, tuple[float, float]] = {
    'D': (8.0, -2.0),   # aspartate — strong repeller within 8 Å
    'E': (9.0, -2.0),   # glutamate
}

# Distance-dependent multipliers for binder residues:
#   ≤4.0 Å  → 1.5×  (direct contact — rare but high-energy)
#   ≤6.0 Å  → 1.2×  (short H-bond range)
#   ≤cutoff → 1.0×  (extended reach / sidechain tip)
_MULT_CLOSE  = 1.5
_MULT_MEDIUM = 1.2
_MULT_FAR    = 1.0


# ---------------------------------------------------------------------------
# PDB parsing
# ---------------------------------------------------------------------------
def parse_pdb(path: str) -> tuple[Optional[tuple[float, float, float]], dict]:
    """Parse a PDB file and return (phosphate_centroid, residue_atoms).

    Parameters
    ----------
    path : str
        Path to PDB (may contain a phospho-Tyr at position 30, with P/OP
        atoms, OR a P4X HETATM ligand split out by split_phospho_to_ligand).

    Returns
    -------
    phos_centroid : (x, y, z) or None
        Centroid of all P/OP atoms found.  None if no phosphate present.
    residues : dict[int, dict[str, tuple]]
        {resnum: {'CA': (x,y,z), 'CB': (x,y,z)}} — only CA and CB stored.
        Hydrogens and digits-prefixed atom names are skipped.
    """
    phos_atoms: list[tuple[float, float, float]] = []
    residues: dict[int, dict] = {}

    with open(path) as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            atom = line[12:16].strip()
            try:
                resnum = int(line[22:26].strip())
                xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            except (ValueError, IndexError):
                continue

            if atom in ('P', 'P1', 'O1P', 'O2P', 'O3P', 'OP1', 'OP2', 'OP3'):
                phos_atoms.append(xyz)
                continue
            if atom.startswith('H') or (atom and atom[0].isdigit()):
                continue
            if atom in ('CA', 'CB'):
                residues.setdefault(resnum, {})[atom] = xyz

    if phos_atoms:
        n = len(phos_atoms)
        centroid: Optional[tuple[float, float, float]] = (
            sum(c[0] for c in phos_atoms) / n,
            sum(c[1] for c in phos_atoms) / n,
            sum(c[2] for c in phos_atoms) / n,
        )
    else:
        centroid = None

    return centroid, residues


def dist(a: tuple, b: tuple) -> float:
    """Euclidean distance between two 3-vectors."""
    return ((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2) ** 0.5


# ---------------------------------------------------------------------------
# Per-residue scoring
# ---------------------------------------------------------------------------
def aa_score(aa: str, d: float) -> tuple[float, str]:
    """Score a single residue at distance *d* from the phosphate centroid.

    Returns
    -------
    score : float
    category : str
        'donor'    — K or R in electrostatic contact
        'h_bond'   — H, S, T, Y, N, Q in H-bond range
        'repeller' — D or E within repulsion cutoff
        'none'     — outside all cutoffs or neutral AA
    """
    if aa in PHOS_BINDERS:
        cutoff, weight = PHOS_BINDERS[aa]
        if d <= cutoff:
            if d <= 4.0:
                mult = _MULT_CLOSE
            elif d <= 6.0:
                mult = _MULT_MEDIUM
            else:
                mult = _MULT_FAR
            category = 'donor' if aa in 'KR' else 'h_bond'
            return weight * mult, category
    elif aa in PHOS_REPELLERS:
        cutoff, weight = PHOS_REPELLERS[aa]
        if d <= cutoff:
            return weight, 'repeller'
    return 0.0, 'none'


# ---------------------------------------------------------------------------
# Full site analysis
# ---------------------------------------------------------------------------
def count_interactions(
    seq: str,
    p_centroid: Optional[tuple],
    residues: dict,
    phos_pos: int = 30,
) -> dict:
    """Analyse all residues in *seq* against a phosphate centroid.

    Excludes the phospho-residue itself (position *phos_pos*).  Tracks
    contacts in three structural regions:

        phos_region (22-29)  — helix immediately N-terminal to phos-Tyr
        central     (31-43)  — central helix / loop
        tail        (44-58)  — hairpin contact face (C-terminal)

    Parameters
    ----------
    seq : str
        59-aa sequence.
    p_centroid : (x, y, z) or None
        Phosphate centroid from parse_pdb().
    residues : dict
        Residue atoms dict from parse_pdb().
    phos_pos : int
        1-indexed position of the phospho-Tyr (default 30).

    Returns
    -------
    dict with keys:
        score, donors, h_bonds, repellers,
        phos_region_donors, central_donors, tail_donors,
        donor_residues (list of "POS AA(dist)" strings),
        repeller_residues (list)
    """
    if p_centroid is None:
        return {
            'score': 0.0,
            'donors': 0, 'h_bonds': 0, 'repellers': 0,
            'donor_residues': [], 'repeller_residues': [],
            'phos_region_donors': 0,
            'central_donors': 0,
            'tail_donors': 0,
        }

    score = 0.0
    donors = h_bonds = repellers = 0
    phos_region = central = tail = 0
    donor_residues: list[str] = []
    repeller_residues: list[str] = []

    for pos, atoms in residues.items():
        if pos < 1 or pos > len(seq) or pos == phos_pos:
            continue
        seq_aa = seq[pos - 1]
        ref = atoms.get('CB') or atoms.get('CA')
        if ref is None:
            continue

        d = dist(ref, p_centroid)
        s, cat = aa_score(seq_aa, d)
        if cat == 'none':
            continue

        score += s

        if cat == 'donor':
            donors += 1
            donor_residues.append(f"{pos}{seq_aa}({d:.1f})")
            if 22 <= pos <= 29:
                phos_region += 1
            elif 31 <= pos <= 43:
                central += 1
            elif 44 <= pos <= 58:
                tail += 1
        elif cat == 'h_bond':
            h_bonds += 1
            if 22 <= pos <= 29:
                phos_region += 1
            elif 31 <= pos <= 43:
                central += 1
            elif 44 <= pos <= 58:
                tail += 1
        elif cat == 'repeller':
            repellers += 1
            repeller_residues.append(f"{pos}{seq_aa}({d:.1f})")

    return {
        'score': round(score, 2),
        'donors': donors,
        'h_bonds': h_bonds,
        'repellers': repellers,
        'donor_residues': donor_residues,
        'repeller_residues': repeller_residues,
        'phos_region_donors': phos_region,
        'central_donors': central,
        'tail_donors': tail,
    }


# ---------------------------------------------------------------------------
# Bidirectional switch score (from mech_iter_bidirectional.py)
# ---------------------------------------------------------------------------
def bidirectional_score(
    seq: str,
    parsed_hairpin: tuple,
    parsed_straight: tuple,
) -> dict:
    """Compute switch metrics simultaneously in both H1 and H2 directions.

    This function does NOT pre-commit to either hypothesis.  It measures how
    much the phosphate-binding differential favours one backbone over the
    other, and returns the BIGGER direction as the winner.

    The effective_score penalises candidates with zero ionic donors (K/R)
    in either backbone (those rely entirely on weak H-bonds) and penalises
    repeller residues that weaken the signal.

    Parameters
    ----------
    seq : str
        59-aa candidate sequence.
    parsed_hairpin : (centroid, residues)
        Output of parse_pdb() for the hairpin backbone (phos_HP).
    parsed_straight : (centroid, residues)
        Output of parse_pdb() for the straight backbone (phos_ST).

    Returns
    -------
    dict with keys:
        switch_magnitude, direction, effective_score,
        hairpin_score, hairpin_donors, hairpin_h_bonds, hairpin_repellers,
        hairpin_tail_contacts, straight_score, straight_donors,
        straight_h_bonds, straight_repellers, straight_phos_contacts,
        hairpin_contacts_str, straight_contacts_str
    """
    hp_p, hp_r = parsed_hairpin
    st_p, st_r = parsed_straight

    hp = count_interactions(seq, hp_p, hp_r)
    st = count_interactions(seq, st_p, st_r)

    diff = hp['score'] - st['score']

    if diff > 0:
        direction = "favors_HAIRPIN"
    elif diff < 0:
        direction = "favors_STRAIGHT"
    else:
        direction = "neutral"

    abs_mag = abs(diff)
    eff_score = abs_mag

    # Penalise designs with no ionic donors (K/R) — they rely on weak H-bonds
    if max(hp['donors'], st['donors']) == 0:
        eff_score -= 3.0

    # Penalise repellers on either backbone
    eff_score -= 0.3 * (hp['repellers'] + st['repellers'])

    hp_contacts = hp.get('donor_residues', []) + hp.get('repeller_residues', [])
    st_contacts = st.get('donor_residues', []) + st.get('repeller_residues', [])

    return {
        'switch_magnitude': round(abs_mag, 2),
        'direction': direction,
        'effective_score': round(eff_score, 2),
        'hairpin_score': hp['score'],
        'hairpin_donors': hp['donors'],
        'hairpin_h_bonds': hp['h_bonds'],
        'hairpin_repellers': hp['repellers'],
        'hairpin_tail_contacts': hp['tail_donors'],
        'straight_score': st['score'],
        'straight_donors': st['donors'],
        'straight_h_bonds': st['h_bonds'],
        'straight_repellers': st['repellers'],
        'straight_phos_contacts': st['phos_region_donors'],
        'hairpin_contacts_str': ' '.join(hp_contacts),
        'straight_contacts_str': ' '.join(st_contacts),
    }


# ---------------------------------------------------------------------------
# CSV-ready mechanism score (from 02_filter_and_mechanism_score.py)
# ---------------------------------------------------------------------------
def mechanism_score(
    seq: str,
    parsed_HP: tuple,
    parsed_ST: tuple,
) -> dict:
    """Compute bidirectional mechanism score and return a flat dict for CSV.

    The direction classification uses donor count differential (integer),
    which is more robust than score differential (float) for ranking.

    Parameters
    ----------
    seq : str
        59-aa candidate.
    parsed_HP : (centroid, residues)
        Hairpin backbone — phos_HP state.
    parsed_ST : (centroid, residues)
        Straight backbone — phos_ST state.

    Returns
    -------
    Flat dict with keys: mech_HP_*, mech_ST_*, mech_donor_diff_HP_minus_ST,
    mech_score_diff_HP_minus_ST, mech_direction.
    """
    hp = count_interactions(seq, parsed_HP[0], parsed_HP[1])
    st = count_interactions(seq, parsed_ST[0], parsed_ST[1])

    if hp['donors'] > st['donors']:
        direction = 'favors_HAIRPIN'
    elif st['donors'] > hp['donors']:
        direction = 'favors_STRAIGHT'
    else:
        direction = 'neutral'

    return {
        'mech_HP_score': hp['score'],
        'mech_HP_donors': hp['donors'],
        'mech_HP_hbonds': hp['h_bonds'],
        'mech_HP_repellers': hp['repellers'],
        'mech_HP_phos_region_donors': hp['phos_region_donors'],
        'mech_HP_central_donors': hp['central_donors'],
        'mech_HP_tail_donors': hp['tail_donors'],
        'mech_HP_donor_str': ' '.join(hp['donor_residues']),

        'mech_ST_score': st['score'],
        'mech_ST_donors': st['donors'],
        'mech_ST_hbonds': st['h_bonds'],
        'mech_ST_repellers': st['repellers'],
        'mech_ST_phos_region_donors': st['phos_region_donors'],
        'mech_ST_central_donors': st['central_donors'],
        'mech_ST_tail_donors': st['tail_donors'],
        'mech_ST_donor_str': ' '.join(st['donor_residues']),

        'mech_donor_diff_HP_minus_ST': hp['donors'] - st['donors'],
        'mech_score_diff_HP_minus_ST': round(hp['score'] - st['score'], 2),
        'mech_direction': direction,
    }


# ---------------------------------------------------------------------------
# H1/H2 classification
# ---------------------------------------------------------------------------
def classify_h1_h2(
    hypothesis: str,
    hp_donors_changed: int,
    st_donors_changed: int,
) -> bool:
    """Return True if the candidate's phosphate contacts match its design hypothesis.

    H1 (phospho stabilises STRAIGHT):
        Hairpin contacts must be BROKEN vs WT.
        Signal: HP donors decreased  OR  ST donors increased.

    H2 (phospho stabilises HAIRPIN):
        Hairpin contacts must be CREATED vs WT.
        Signal: HP donors increased  OR  ST donors decreased.

    Parameters
    ----------
    hypothesis : str
        'H1' or 'H2'.
    hp_donors_changed : int
        mech_HP_donors(design) − mech_HP_donors(WT).
    st_donors_changed : int
        mech_ST_donors(design) − mech_ST_donors(WT).
    """
    if hypothesis == 'H1':
        return (hp_donors_changed < 0) or (st_donors_changed > 0)
    if hypothesis == 'H2':
        return (hp_donors_changed > 0) or (st_donors_changed < 0)
    return False
