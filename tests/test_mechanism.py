"""
test_mechanism.py — unit tests for phosphate contact scoring.

Tests the bidirectional_score function with synthetic toy PDB coordinates
so that no real PDB files are required for CI.
"""

import math
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from phosphoswitch_design.mechanism import (
    PHOS_BINDERS,
    PHOS_REPELLERS,
    dist,
    aa_score,
    count_interactions,
    bidirectional_score,
    mechanism_score,
    classify_h1_h2,
)


# ---------------------------------------------------------------------------
# Toy geometry: place phosphate centroid at origin, residues at fixed distances
# ---------------------------------------------------------------------------
PHOS_AT_ORIGIN = (0.0, 0.0, 0.0)

# Residue at 5 Å — within K/R/H/N/Q cutoffs but outside S/T (6 Å)
RES_5A = (5.0, 0.0, 0.0)

# Residue at 3.5 Å — close contact (mult 1.5×)
RES_3A = (3.5, 0.0, 0.0)

# Residue at 12 Å — outside all binder cutoffs but within R (11 Å)? No: 12>11
RES_12A = (12.0, 0.0, 0.0)

# Residue at 7.5 Å — within D/E repeller (8/9 Å) cutoff
RES_7A = (7.5, 0.0, 0.0)


def _make_parsed(centroid, residue_positions: dict[int, tuple]) -> tuple:
    """Build a (centroid, residues) tuple for testing count_interactions."""
    residues = {pos: {'CA': xyz} for pos, xyz in residue_positions.items()}
    return centroid, residues


class TestDist(unittest.TestCase):
    def test_origin(self):
        self.assertAlmostEqual(dist((0,0,0), (0,0,0)), 0.0)

    def test_unit(self):
        self.assertAlmostEqual(dist((0,0,0), (1,0,0)), 1.0)

    def test_3d(self):
        d = dist((1,2,3), (4,6,3))
        self.assertAlmostEqual(d, 5.0)


class TestAaScore(unittest.TestCase):
    def test_lysine_medium_range(self):
        # K at 5 Å: within cutoff 10 Å, medium range (≤6 Å), weight 1.5, mult 1.2
        score, cat = aa_score('K', 5.0)
        self.assertEqual(cat, 'donor')
        self.assertAlmostEqual(score, 1.5 * 1.2)

    def test_lysine_close(self):
        # K at 3.5 Å: close contact, mult 1.5
        score, cat = aa_score('K', 3.5)
        self.assertEqual(cat, 'donor')
        self.assertAlmostEqual(score, 1.5 * 1.5)

    def test_arginine_far(self):
        # R at 10.5 Å: within cutoff 11 Å, far range (>6 Å), weight 2.0, mult 1.0
        score, cat = aa_score('R', 10.5)
        self.assertEqual(cat, 'donor')
        self.assertAlmostEqual(score, 2.0 * 1.0)

    def test_arginine_outside_cutoff(self):
        score, cat = aa_score('R', 12.0)
        self.assertEqual(cat, 'none')
        self.assertEqual(score, 0.0)

    def test_histidine_hbond(self):
        score, cat = aa_score('H', 7.0)
        self.assertEqual(cat, 'h_bond')
        self.assertGreater(score, 0)

    def test_serine_too_far(self):
        # S cutoff is 6.0 Å; at 7 Å → none
        score, cat = aa_score('S', 7.0)
        self.assertEqual(cat, 'none')
        self.assertEqual(score, 0.0)

    def test_aspartate_repeller(self):
        # D at 7.5 Å: within cutoff 8 Å → repeller
        score, cat = aa_score('D', 7.5)
        self.assertEqual(cat, 'repeller')
        self.assertLess(score, 0)

    def test_glutamate_outside_cutoff(self):
        # E cutoff 9 Å; at 10 Å → none
        score, cat = aa_score('E', 10.0)
        self.assertEqual(cat, 'none')

    def test_alanine_neutral(self):
        score, cat = aa_score('A', 3.0)
        self.assertEqual(cat, 'none')
        self.assertEqual(score, 0.0)


class TestCountInteractions(unittest.TestCase):
    def test_no_centroid(self):
        """If no phosphate centroid, all counts should be zero."""
        parsed = _make_parsed(None, {25: (3.0, 0.0, 0.0)})
        result = count_interactions("A" * 59, None, parsed[1])
        self.assertEqual(result['score'], 0)
        self.assertEqual(result['donors'], 0)

    def test_single_lysine_donor(self):
        # Position 25, K at 5 Å from phos centroid at origin
        # Sequence: 59 aa, pos 25 = K
        seq = "A" * 24 + "K" + "A" * 34
        parsed = _make_parsed(PHOS_AT_ORIGIN, {25: RES_5A})
        result = count_interactions(seq, PHOS_AT_ORIGIN, parsed[1])
        self.assertEqual(result['donors'], 1)
        self.assertGreater(result['score'], 0)
        # Position 25 is in range 22-29 (phos_region)
        self.assertEqual(result['phos_region_donors'], 1)

    def test_skip_phos_position(self):
        # Position 30 is excluded (the phospho-Tyr itself)
        seq = "A" * 29 + "Y" + "A" * 29
        parsed = _make_parsed(PHOS_AT_ORIGIN, {30: (2.0, 0.0, 0.0)})
        result = count_interactions(seq, PHOS_AT_ORIGIN, parsed[1])
        self.assertEqual(result['donors'], 0)
        self.assertEqual(result['h_bonds'], 0)

    def test_repeller(self):
        seq = "A" * 31 + "D" + "A" * 27
        parsed = _make_parsed(PHOS_AT_ORIGIN, {32: RES_7A})
        result = count_interactions(seq, PHOS_AT_ORIGIN, parsed[1])
        self.assertEqual(result['repellers'], 1)
        self.assertLess(result['score'], 0)

    def test_tail_region(self):
        # Position 50 is in the tail (44-58)
        seq = "A" * 49 + "K" + "A" * 9
        parsed = _make_parsed(PHOS_AT_ORIGIN, {50: RES_5A})
        result = count_interactions(seq, PHOS_AT_ORIGIN, parsed[1])
        self.assertEqual(result['tail_donors'], 1)


class TestBidirectionalScore(unittest.TestCase):
    def _make_hp_st_favoring_hp(self):
        """Make geometry where HP has K close, ST has K far → favors HP."""
        # HP: K at position 50, distance 5 Å
        hp_seq = "A" * 49 + "K" + "A" * 9
        parsed_hp = _make_parsed(PHOS_AT_ORIGIN, {50: RES_5A})

        # ST: K at position 50, distance 12 Å (outside cutoff)
        parsed_st = _make_parsed(PHOS_AT_ORIGIN, {50: RES_12A})
        return hp_seq, parsed_hp, parsed_st

    def test_favors_hairpin(self):
        seq, parsed_hp, parsed_st = self._make_hp_st_favoring_hp()
        result = bidirectional_score(seq, parsed_hp, parsed_st)
        self.assertEqual(result['direction'], 'favors_HAIRPIN')
        self.assertGreater(result['switch_magnitude'], 0)

    def test_favors_straight(self):
        seq = "A" * 49 + "K" + "A" * 9
        # Swap: ST close, HP far
        parsed_hp = _make_parsed(PHOS_AT_ORIGIN, {50: RES_12A})
        parsed_st = _make_parsed(PHOS_AT_ORIGIN, {50: RES_5A})
        result = bidirectional_score(seq, parsed_hp, parsed_st)
        self.assertEqual(result['direction'], 'favors_STRAIGHT')

    def test_neutral_no_contacts(self):
        # Alanine everywhere → no contacts on either backbone
        seq = "A" * 59
        parsed_hp = _make_parsed(PHOS_AT_ORIGIN, {50: RES_5A})
        parsed_st = _make_parsed(PHOS_AT_ORIGIN, {50: RES_5A})
        result = bidirectional_score(seq, parsed_hp, parsed_st)
        self.assertEqual(result['direction'], 'neutral')
        self.assertEqual(result['switch_magnitude'], 0)

    def test_zero_donors_penalty(self):
        # H-bond donors only → effective_score penalised by 3.0
        seq = "A" * 24 + "S" + "A" * 34  # S at pos 25 (within 6 Å)
        parsed_hp = _make_parsed(PHOS_AT_ORIGIN, {25: (4.0, 0.0, 0.0)})
        parsed_st = _make_parsed(PHOS_AT_ORIGIN, {25: RES_12A})
        result = bidirectional_score(seq, parsed_hp, parsed_st)
        # No K or R → donors=0 → -3.0 penalty
        self.assertLess(result['effective_score'], result['switch_magnitude'])

    def test_repeller_penalty(self):
        # D repeller → negative contribution to effective_score
        seq = "A" * 31 + "D" + "A" * 27
        parsed_hp = _make_parsed(PHOS_AT_ORIGIN, {32: RES_7A})
        parsed_st = _make_parsed(PHOS_AT_ORIGIN, {32: RES_7A})
        result = bidirectional_score(seq, parsed_hp, parsed_st)
        self.assertLess(result['effective_score'], 0)


class TestMechanismScore(unittest.TestCase):
    def test_keys_present(self):
        """mechanism_score() must return all expected CSV columns."""
        seq = "A" * 59
        parsed = _make_parsed(PHOS_AT_ORIGIN, {})
        result = mechanism_score(seq, parsed, parsed)
        required_keys = [
            'mech_HP_score', 'mech_HP_donors', 'mech_HP_hbonds',
            'mech_ST_score', 'mech_ST_donors', 'mech_ST_hbonds',
            'mech_donor_diff_HP_minus_ST', 'mech_score_diff_HP_minus_ST',
            'mech_direction',
        ]
        for k in required_keys:
            self.assertIn(k, result, f"Missing key: {k}")

    def test_neutral_when_symmetric(self):
        """Same contacts on both backbones → neutral direction."""
        seq = "A" * 49 + "K" + "A" * 9
        parsed = _make_parsed(PHOS_AT_ORIGIN, {50: RES_5A})
        result = mechanism_score(seq, parsed, parsed)
        self.assertEqual(result['mech_direction'], 'neutral')
        self.assertEqual(result['mech_donor_diff_HP_minus_ST'], 0)


class TestClassifyH1H2(unittest.TestCase):
    def test_h1_hp_decreased(self):
        # H1: HP donors decreased → True
        self.assertTrue(classify_h1_h2('H1', hp_donors_changed=-1, st_donors_changed=0))

    def test_h1_st_increased(self):
        # H1: ST donors increased → True
        self.assertTrue(classify_h1_h2('H1', hp_donors_changed=0, st_donors_changed=1))

    def test_h1_both_zero(self):
        # H1: no change → False
        self.assertFalse(classify_h1_h2('H1', hp_donors_changed=0, st_donors_changed=0))

    def test_h2_hp_increased(self):
        # H2: HP donors increased → True
        self.assertTrue(classify_h1_h2('H2', hp_donors_changed=1, st_donors_changed=0))

    def test_h2_st_decreased(self):
        # H2: ST donors decreased → True
        self.assertTrue(classify_h1_h2('H2', hp_donors_changed=0, st_donors_changed=-1))

    def test_h2_wrong_direction(self):
        # H2: HP donors decreased (wrong for H2) → False if ST unchanged too
        self.assertFalse(classify_h1_h2('H2', hp_donors_changed=-1, st_donors_changed=0))

    def test_unknown_hypothesis(self):
        self.assertFalse(classify_h1_h2('H3', hp_donors_changed=1, st_donors_changed=-1))


if __name__ == "__main__":
    unittest.main(verbosity=2)
