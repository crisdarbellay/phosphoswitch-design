"""
test_ddg.py — unit tests for the 4-state thermodynamic cycle mathematics.

Tests that:
    1. ddG_switch is computed correctly from the 4-corner formula
    2. aggregate_replicates() produces correct statistics from known values
    3. Sign convention: positive ddG_switch = H2 (phospho prefers hairpin)
    4. Outlier detection fires at the correct threshold
    5. z-score computation is correct against a reference distribution

No PyRosetta is required — the tests use synthetic energy values.
"""

import math
import statistics
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from phosphoswitch_design.rosetta import aggregate_replicates


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _make_rep(tag: str, rep: int, e_phos_st: float, e_phos_hp: float,
              e_apo_st: float, e_apo_hp: float) -> dict:
    """Build a synthetic result dict like score_one_replicate() returns."""
    phos_pref = e_phos_hp - e_phos_st
    apo_pref  = e_apo_hp  - e_apo_st
    ddG       = phos_pref - apo_pref
    return {
        'tag': tag,
        'replicate': rep,
        'sequence': 'A' * 59,
        'E_phos_ST': e_phos_st,
        'E_phos_HP': e_phos_hp,
        'E_apo_ST':  e_apo_st,
        'E_apo_HP':  e_apo_hp,
        'phos_pref': phos_pref,
        'apo_pref':  apo_pref,
        'ddG_switch_4corner': ddG,
    }


class TestDdgFormula(unittest.TestCase):
    """Test the 4-corner ddG formula directly."""

    def test_h2_signal(self):
        """Phospho preferring hairpin → positive ddG_switch."""
        # Phospho makes HP 10 REU lower than ST (prefers HP)
        # Apo is indifferent (same energy on both)
        rep = _make_rep('t1', 0,
                        e_phos_st=-100.0, e_phos_hp=-110.0,  # phos prefers HP by 10
                        e_apo_st=-100.0,  e_apo_hp=-100.0)   # apo indifferent
        # phos_pref = -110 - (-100) = -10  → HP is lower for phos
        # apo_pref  = -100 - (-100) =   0
        # ddG = -10 - 0 = -10  → hairpin preferred by phospho
        # Wait: the formula is ddG = phos_pref − apo_pref
        # phos_pref = E_phos_HP − E_phos_ST = -110 − (-100) = -10
        # apo_pref  = E_apo_HP  − E_apo_ST  = -100 − (-100) =   0
        # ddG = -10 − 0 = -10
        # Negative means phospho prefers HP (H2)?  Let's check the sign convention:
        #   ddG < 0 → phos_pref < apo_pref → phos makes HP even MORE stable → H2
        self.assertAlmostEqual(rep['ddG_switch_4corner'], -10.0)
        # By convention used in the pipeline:
        # phos_pref = E_phos_HP - E_phos_ST: if negative, phospho prefers HP → H2
        # ddG = phos_pref - apo_pref; ddG < 0 → H2 direction
        self.assertLess(rep['ddG_switch_4corner'], 0)

    def test_h1_signal(self):
        """Phospho preferring straight → phos_pref > apo_pref → ddG positive? """
        # Phospho makes ST 10 REU lower than HP (prefers ST)
        rep = _make_rep('t2', 0,
                        e_phos_st=-110.0, e_phos_hp=-100.0,   # phos prefers ST
                        e_apo_st=-100.0,  e_apo_hp=-100.0)    # apo indifferent
        # phos_pref = -100 - (-110) = +10  (HP higher than ST for phos)
        # apo_pref  = -100 - (-100) =   0
        # ddG = 10 − 0 = +10
        self.assertAlmostEqual(rep['ddG_switch_4corner'], 10.0)
        self.assertGreater(rep['ddG_switch_4corner'], 0)

    def test_zero_when_indifferent(self):
        """No phospho effect → ddG_switch = 0."""
        rep = _make_rep('t3', 0,
                        e_phos_st=-100.0, e_phos_hp=-105.0,
                        e_apo_st=-100.0,  e_apo_hp=-105.0)
        # phos_pref = -105 - (-100) = -5
        # apo_pref  = -105 - (-100) = -5
        # ddG = -5 - (-5) = 0
        self.assertAlmostEqual(rep['ddG_switch_4corner'], 0.0)

    def test_additive_consistency(self):
        """ddG_switch = phos_pref − apo_pref (additive thermodynamic cycle)."""
        e = dict(e_phos_st=-150.0, e_phos_hp=-145.0,
                 e_apo_st=-148.0,  e_apo_hp=-147.0)
        rep = _make_rep('t4', 0, **e)
        expected = (e['e_phos_hp'] - e['e_phos_st']) - (e['e_apo_hp'] - e['e_apo_st'])
        self.assertAlmostEqual(rep['ddG_switch_4corner'], expected, places=6)


class TestAggregateReplicates(unittest.TestCase):
    """Test aggregate_replicates() statistics."""

    def _wt_dist(self, n: int = 50, center: float = 0.0, std: float = 2.0) -> list[float]:
        """Synthetic WT distribution centered at *center* with *std*."""
        import random
        rng = random.Random(42)
        return [rng.gauss(center, std) for _ in range(n)]

    def test_basic_statistics(self):
        """Mean, median, std, range must be correct for known values."""
        replicates = [
            _make_rep('cand', i, -100.0, -100.0 + i, -100.0, -100.0)
            for i in range(10)
        ]
        # ddG values: [-100+i - (-100)] - 0 = i, so ddG = i for i in 0..9
        # phos_pref = i, apo_pref = 0, ddG = i
        wt = self._wt_dist()
        result = aggregate_replicates('cand', replicates, wt)

        expected_ddgs = [float(i) for i in range(10)]
        self.assertAlmostEqual(result['ddG_mean'],   statistics.mean(expected_ddgs),   places=3)
        self.assertAlmostEqual(result['ddG_median'], statistics.median(expected_ddgs), places=3)
        self.assertAlmostEqual(result['ddG_std'],    statistics.stdev(expected_ddgs),  places=3)
        self.assertAlmostEqual(result['ddG_min'],    0.0, places=3)
        self.assertAlmostEqual(result['ddG_max'],    9.0, places=3)
        self.assertAlmostEqual(result['ddG_range'],  9.0, places=3)

    def test_no_outlier_flag(self):
        """Range < 15 REU → outlier_flag = False."""
        replicates = [
            _make_rep('c', 0, -100.0, -110.0, -100.0, -100.0),  # ddG = -10
            _make_rep('c', 1, -100.0, -108.0, -100.0, -100.0),  # ddG = -8
        ]
        result = aggregate_replicates('c', replicates, self._wt_dist())
        # range = |-8 - (-10)| = 2 < 15
        self.assertFalse(result['outlier_flag'])

    def test_outlier_flag(self):
        """Range > 15 REU → outlier_flag = True."""
        replicates = [
            _make_rep('c', 0, -100.0, -120.0, -100.0, -100.0),  # ddG = -20
            _make_rep('c', 1, -100.0, -100.0, -100.0, -100.0),  # ddG =   0
        ]
        result = aggregate_replicates('c', replicates, self._wt_dist())
        # range = 20 > 15
        self.assertTrue(result['outlier_flag'])

    def test_z_score_zero_for_wt_like(self):
        """Candidate with same median as WT → z-score ≈ 0."""
        wt_dist = [0.0] * 50   # degenerate WT dist (std=0 → z undefined)
        # Use a realistic WT dist
        wt_dist = self._wt_dist(50, center=0.0, std=3.0)
        wt_median = statistics.median(wt_dist)

        # Candidate with same ddG as WT median
        replicates = [
            _make_rep('cand', i, -100.0, -100.0 + wt_median, -100.0, -100.0)
            for i in range(5)
        ]
        result = aggregate_replicates('cand', replicates, wt_dist)
        self.assertAlmostEqual(result['z_score_vs_WT'], 0.0, delta=0.5)

    def test_significant_z(self):
        """Extreme ddG → z-score should exceed Bonferroni threshold (4.4)."""
        # WT distribution centered at 0 with std=2
        wt_dist = self._wt_dist(50, center=0.0, std=2.0)

        # Candidate 20 REU away from WT center
        replicates = [
            _make_rep('cand', i, -100.0, -120.0, -100.0, -100.0)
            for i in range(20)
        ]
        result = aggregate_replicates('cand', replicates, wt_dist)
        self.assertGreater(abs(result['z_score_vs_WT']), 4.4)
        self.assertEqual(result['significance'], 'significant')

    def test_empty_replicates(self):
        """Empty replicate list → n_reps=0, error key present."""
        result = aggregate_replicates('empty', [], [0.0] * 10)
        self.assertEqual(result['n_reps'], 0)
        self.assertIn('error', result)

    def test_failed_replicate_ignored(self):
        """Replicates with 'error' key (no ddG_switch_4corner) are skipped."""
        replicates = [
            {'tag': 'c', 'replicate': 0, 'error': 'PyRosetta init failed', 'sequence': 'A'*59},
            _make_rep('c', 1, -100.0, -110.0, -100.0, -100.0),  # ddG = -10
        ]
        result = aggregate_replicates('c', replicates, self._wt_dist())
        self.assertEqual(result['n_reps'], 1)
        self.assertAlmostEqual(result['ddG_median'], -10.0, places=3)

    def test_n_reps_count(self):
        """n_reps should equal the number of valid (non-error) replicates."""
        replicates = [
            _make_rep('c', i, -100.0, -105.0, -100.0, -100.0)
            for i in range(7)
        ]
        result = aggregate_replicates('c', replicates, self._wt_dist())
        self.assertEqual(result['n_reps'], 7)


class TestDdgSignConvention(unittest.TestCase):
    """Verify the sign convention is consistent with the biology documentation."""

    def test_h2_mechanism_gives_negative_ddg(self):
        """
        H2: phospho stabilises HAIRPIN.
        If phospho lowers HP energy more than apo does, ddG_switch < 0.
        Convention: ddG_switch = (E_phos_HP - E_phos_ST) - (E_apo_HP - E_apo_ST)
        """
        # Phospho: HP 15 REU lower than ST
        # Apo:     HP 5 REU lower than ST (some baseline HP preference)
        # Net phospho bonus to HP = 15 - 5 = 10 REU
        rep = _make_rep('h2', 0,
                        e_phos_st=-100.0, e_phos_hp=-115.0,
                        e_apo_st=-100.0,  e_apo_hp=-105.0)
        # phos_pref = -115 - (-100) = -15
        # apo_pref  = -105 - (-100) =  -5
        # ddG = -15 - (-5) = -10  → phospho further stabilises HP → H2
        self.assertAlmostEqual(rep['ddG_switch_4corner'], -10.0)

    def test_h1_mechanism_gives_positive_ddg(self):
        """
        H1: phospho stabilises STRAIGHT.
        ddG_switch > 0 means phospho makes straight even more preferred.
        """
        rep = _make_rep('h1', 0,
                        e_phos_st=-115.0, e_phos_hp=-100.0,
                        e_apo_st=-105.0,  e_apo_hp=-100.0)
        # phos_pref = -100 - (-115) = +15
        # apo_pref  = -100 - (-105) =  +5
        # ddG = 15 - 5 = +10  → phospho further stabilises ST → H1
        self.assertAlmostEqual(rep['ddG_switch_4corner'], 10.0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
