#!/usr/bin/env python3
"""
03_select.py — top diverse candidate selection for Rosetta validation.

Selects ~10,000 candidates from the ~100k plausible+scored pool using:

    1. H1/H2 directional contact filter
       Eliminates designs whose phosphate-contact profile contradicts the
       design hypothesis.  This is the most important gate before expensive
       Rosetta validation.

    2. Rank by |deviation_score_diff| within each direction pool.
       Largest deviation from WT signal → strongest switch candidates.

    3. Diversity caps
       Max N per subspace, max N per mutation signature.
       Direction allocation: 40% HP-favoring, 40% ST-favoring, 20% neutral.

Output: ~10,000 high-quality candidates ready for deep Rosetta validation
(stage 04 runs N=20 FastRelax replicates per candidate, ~25h with 16 workers).

Estimated Rosetta compute time at N=10,000, N=20 reps:
    8  workers: ~83 h
    16 workers: ~42 h
    24 workers: ~28 h

Usage
-----
    python scripts/03_select.py --help
    python scripts/03_select.py --input-csv outputs/02_plausible_with_mechanism.csv
    python scripts/03_select.py --n-top 5000 --no-hypothesis-filter
"""

import argparse
import csv
import os
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from phosphoswitch_design.filtering import (
    passes_h1_h2_filter,
    select_with_diversity,
)
from phosphoswitch_design.io_utils import write_csv


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--input-csv",  default="outputs/02_plausible_with_mechanism.csv")
    ap.add_argument("--out-csv",    default="outputs/03_top_diverse_candidates.csv")
    ap.add_argument("--n-top",      type=int, default=10000,
                    help="Total candidates to select (default 10,000)")
    ap.add_argument("--max-per-subspace",   type=int, default=1000)
    ap.add_argument("--max-per-signature",  type=int, default=30)
    ap.add_argument("--frac-hp",    type=float, default=0.4,
                    help="Fraction HP-favoring candidates (default 0.4)")
    ap.add_argument("--frac-st",    type=float, default=0.4,
                    help="Fraction ST-favoring candidates (default 0.4)")
    ap.add_argument("--no-hypothesis-filter", action="store_true",
                    help="Disable H1/H2 directional contact filter")
    return ap


def main() -> None:
    ap = build_parser()
    args = ap.parse_args()

    if not os.path.exists(args.input_csv):
        sys.exit(f"ERROR: {args.input_csv} not found.\nRun 02_filter.py first.")

    print("=" * 78)
    print("  03_select — diverse top candidate selection")
    print("=" * 78)
    if args.no_hypothesis_filter:
        print("  NOTE: H1/H2 directional contact filter DISABLED")

    # ------------------------------------------------------------------
    # Load CSV
    # ------------------------------------------------------------------
    all_records = []
    template_row = None

    with open(args.input_csv) as fh:
        for row in csv.DictReader(fh):
            try:
                row['_dev']          = float(row.get('deviation_score_diff', 0))
                row['_n_muts']       = int(row.get('n_mutations', 99))
                row['_hp_repellers'] = int(row.get('mech_HP_repellers', 99))
                row['_hyp_match']    = passes_h1_h2_filter(row)
            except (ValueError, KeyError):
                continue

            if row.get('run') == 'TEMPLATE':
                template_row = row
            else:
                all_records.append(row)

    print(f"\n  Loaded {len(all_records):,} candidates + 1 WT template")

    # ------------------------------------------------------------------
    # H1/H2 directional contact filter
    # ------------------------------------------------------------------
    if not args.no_hypothesis_filter:
        before = len(all_records)
        all_records = [r for r in all_records if r['_hyp_match']]
        after = len(all_records)
        pct = 100 * after / before if before else 0
        print(f"\n  H1/H2 filter: {before:,} → {after:,} ({pct:.1f}% pass)")
        if after == 0:
            print("  WARNING: no candidates passed the filter.")
            print(f"  Check {args.input_csv} for non-WT rows.")
            os.makedirs(os.path.dirname(os.path.abspath(args.out_csv)), exist_ok=True)
            with open(args.out_csv, 'w') as fh:
                fh.write("tag,sequence,run,subspace,hypothesis,_selection_pool\n")
            return

    # ------------------------------------------------------------------
    # Direction split by deviation_score_diff
    # ------------------------------------------------------------------
    favor_HP = [r for r in all_records if r['_dev'] > 1.0]
    favor_ST = [r for r in all_records if r['_dev'] < -1.0]
    neutral  = [r for r in all_records if abs(r['_dev']) <= 1.0]

    favor_HP.sort(key=lambda r: -r['_dev'])
    favor_ST.sort(key=lambda r:  r['_dev'])
    neutral.sort( key=lambda r: -abs(r['_dev']))

    print(f"\n  Direction split:")
    print(f"    favor_HAIRPIN  (dev >+1): {len(favor_HP):>8,}")
    print(f"    favor_STRAIGHT (dev <-1): {len(favor_ST):>8,}")
    print(f"    neutral        (|dev|≤1): {len(neutral):>8,}")

    # ------------------------------------------------------------------
    # Diversity-capped selection
    # ------------------------------------------------------------------
    n_hp  = int(args.n_top * args.frac_hp)
    n_st  = int(args.n_top * args.frac_st)
    n_neu = args.n_top - n_hp - n_st

    print(f"\n  Selecting with caps: max_per_subspace={args.max_per_subspace}, "
          f"max_per_signature={args.max_per_signature}")

    hp_sel, _, _  = select_with_diversity(favor_HP, n_hp,  'favor_HP',
                                           args.max_per_signature, args.max_per_subspace)
    st_sel, _, _  = select_with_diversity(favor_ST, n_st,  'favor_ST',
                                           args.max_per_signature, args.max_per_subspace)
    neu_sel, _, _ = select_with_diversity(neutral,  n_neu, 'neutral',
                                           args.max_per_signature, args.max_per_subspace)

    print(f"    HP-favoring: {len(hp_sel):>6,} / {n_hp}")
    print(f"    ST-favoring: {len(st_sel):>6,} / {n_st}")
    print(f"    neutral:     {len(neu_sel):>6,} / {n_neu}")

    selected = hp_sel + st_sel + neu_sel
    if template_row:
        template_row['_selection_pool'] = 'WT_baseline'
        selected = [template_row] + selected

    print(f"\n  Total selected: {len(selected):,} ({len(selected)-1:,} + 1 WT)")

    # Selection breakdown
    by_track    = defaultdict(int)
    by_subspace = defaultdict(int)
    by_hyp      = defaultdict(int)
    for r in selected:
        if r.get('run') == 'TEMPLATE':
            continue
        by_track[r.get('run', '?')] += 1
        by_subspace[r.get('subspace', '?')] += 1
        by_hyp[r.get('hypothesis', '?')] += 1

    print(f"\n  By track:      { dict(sorted(by_track.items())) }")
    print(f"  By hypothesis: { dict(sorted(by_hyp.items())) }")

    # Rosetta time estimate
    n = len(selected)
    print(f"\n  Estimated Rosetta time (N=20 reps × 4 corners, ~2s each):")
    for w in [8, 16, 24]:
        t = n * 20 * 4 * 2 / w / 60
        print(f"    {w} workers: {t:.0f} min ({t/60:.1f}h)")

    # ------------------------------------------------------------------
    # Write output
    # ------------------------------------------------------------------
    if selected:
        all_keys: set[str] = set()
        for r in selected:
            all_keys.update(r.keys())
        cols = sorted(c for c in all_keys if not c.startswith('_'))
        cols.append('_selection_pool')
        write_csv(args.out_csv, selected, fieldnames=cols)
        print(f"\n  Output: {args.out_csv}")

    print(f"\n  Next step:")
    print(f"    python scripts/04_rosetta.py --input-csv {args.out_csv}")


if __name__ == "__main__":
    main()
