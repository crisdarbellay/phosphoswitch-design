#!/usr/bin/env python3
"""
04_rosetta.py — deep PyRosetta 4-state thermodynamic validation.

Scores each candidate on the full 4-corner thermodynamic cycle:

    phos_ST (phospho + straight)    phos_HP (phospho + hairpin)
         \\                              /
          \\                            /
    apo_ST  (apo + straight)      apo_HP  (apo + hairpin)

    ddG_switch = (E_phos_HP - E_phos_ST) - (E_apo_HP - E_apo_ST)

    ddG_switch > 0  →  phospho increases HAIRPIN preference  (H2 mechanism)
    ddG_switch < 0  →  phospho increases STRAIGHT preference (H1 mechanism)
    ddG_switch ≈ 0  →  phospho has no conformational effect

Statistics (N=20 replicates per candidate)
------------------------------------------
Each replicate uses a different stochastic seed derived from (replicate, tag).
FastRelax: 3 cycles, backbone constrained, sidechains free.

From the replicate distribution we report:
    mean, median, std, range, z-score vs WT, Bonferroni significance,
    outlier_flag (range > 15 REU = likely stochastic explosion).

WT is scored with N=50 replicates to establish a tight reference distribution.

Known artifacts
---------------
- Histidine at fixed protonation state inflates score artificially
  (real His mostly neutral at pH 7-7.4)
- FastRelax noise floor: ~±3-5 REU per replicate → ±7 REU for 4-corner cycle
  → no single-shot result below ~10 REU should be trusted

Resume
------
Per-replicate JSON cache files are written to --work-dir.  Re-running with
--resume skips any (tag, replicate) with a valid cached JSON.

Usage
-----
    python scripts/04_rosetta.py --help
    python scripts/04_rosetta.py --workers 16
    python scripts/04_rosetta.py --n-reps 5 --workers 4  # quick test
    python scripts/04_rosetta.py --resume
"""

import argparse
import csv
import json
import os
import sys
import time
from collections import defaultdict
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from phosphoswitch_design.rosetta import (
    score_one_replicate,
    aggregate_replicates,
    make_apo_pdb,
)
from phosphoswitch_design.io_utils import write_csv


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--input-csv",  default="outputs/03_top_diverse_candidates.csv")
    ap.add_argument("--out-csv",    default="outputs/04_deep_validated.csv")
    ap.add_argument("--work-dir",   default="outputs/04_rosetta_cache",
                    help="Directory for per-replicate JSON cache files")
    ap.add_argument("--phos-st",    default="output/phase1/stateA_phospho.pdb",
                    help="Straight backbone + phosphate PDB (phos_ST)")
    ap.add_argument("--phos-hp",    default="output/phase1/stateA_phospho_pulled.pdb",
                    help="Hairpin backbone + phosphate PDB (phos_HP)")
    ap.add_argument("--apo-st",     default="outputs/04_rosetta_cache/apo_ST.pdb",
                    help="Apo straight PDB (auto-derived if absent)")
    ap.add_argument("--apo-hp",     default="outputs/04_rosetta_cache/apo_HP.pdb",
                    help="Apo hairpin PDB (auto-derived if absent)")
    ap.add_argument("--n-reps",     type=int, default=20,
                    help="Replicates per candidate (default 20)")
    ap.add_argument("--n-reps-wt",  type=int, default=50,
                    help="Replicates for WT baseline (default 50)")
    ap.add_argument("--workers",    type=int, default=8,
                    help="Parallel PyRosetta workers (default 8)")
    ap.add_argument("--outlier-threshold", type=float, default=15.0,
                    help="ddG range threshold for outlier flag (REU, default 15)")
    ap.add_argument("--resume",     action="store_true",
                    help="Skip (tag, rep) pairs with cached JSON")
    ap.add_argument("--max-candidates", type=int, default=None,
                    help="Limit candidates for testing")
    return ap


def load_candidates(csv_path: str, max_candidates: int | None = None) -> list[dict]:
    """Load candidates from stage 03 CSV.  Returns list of {tag, sequence, ...}."""
    candidates = []
    with open(csv_path) as fh:
        for row in csv.DictReader(fh):
            seq = row.get('sequence', '')
            if not seq or len(seq) != 59:
                continue
            tag = (f"{row.get('run','?')}_{row.get('subspace','?')}_"
                   f"T{row.get('temperature','?')}_{len(candidates)}")
            row['tag'] = tag
            candidates.append(row)
            if max_candidates and len(candidates) >= max_candidates:
                break
    return candidates


def main() -> None:
    ap = build_parser()
    args = ap.parse_args()

    if not os.path.exists(args.input_csv):
        sys.exit(f"ERROR: {args.input_csv} not found. Run 03_select.py first.")

    for pdb_path, label in [(args.phos_st, "--phos-st"), (args.phos_hp, "--phos-hp")]:
        if not os.path.exists(pdb_path):
            sys.exit(f"ERROR: backbone PDB not found: {pdb_path}\n"
                     f"Set {label} to your phase1 PDB path.")

    os.makedirs(args.work_dir, exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(args.out_csv)), exist_ok=True)

    # Derive apo PDBs
    make_apo_pdb(args.phos_st, args.apo_st)
    make_apo_pdb(args.phos_hp, args.apo_hp)

    pdbs = {
        'phos_ST': args.phos_st,
        'phos_HP': args.phos_hp,
        'apo_ST':  args.apo_st,
        'apo_HP':  args.apo_hp,
    }

    print("=" * 78)
    print("  04_rosetta — deep 4-state thermodynamic validation")
    print("=" * 78)
    print(f"  phos_ST: {args.phos_st}")
    print(f"  phos_HP: {args.phos_hp}")
    print(f"  apo_ST:  {args.apo_st}  (auto-derived)")
    print(f"  apo_HP:  {args.apo_hp}  (auto-derived)")
    print(f"  Replicates per candidate: {args.n_reps}")
    print(f"  Replicates for WT:        {args.n_reps_wt}")
    print(f"  Workers:                  {args.workers}")
    print(f"  Resume:                   {args.resume}")

    # ------------------------------------------------------------------
    # Load candidates
    # ------------------------------------------------------------------
    candidates = load_candidates(args.input_csv, args.max_candidates)
    wt_candidates = [c for c in candidates if c.get('run') == 'TEMPLATE']
    design_candidates = [c for c in candidates if c.get('run') != 'TEMPLATE']

    n_reps_wt = args.n_reps_wt if wt_candidates else 0
    n_combos = len(design_candidates) * args.n_reps + len(wt_candidates) * n_reps_wt
    est_h = n_combos * 2 / args.workers / 3600

    print(f"\n  Designs:    {len(design_candidates):,}")
    print(f"  WT:         {len(wt_candidates)}")
    print(f"  Total jobs: {n_combos:,}")
    print(f"  Est time:   ~{est_h:.1f}h")

    # ------------------------------------------------------------------
    # Build job list
    # ------------------------------------------------------------------
    def make_jobs(cands, n_reps):
        return [
            (c, rep, args.work_dir, pdbs)
            for c in cands
            for rep in range(n_reps)
        ]

    jobs = make_jobs(wt_candidates, n_reps_wt) + make_jobs(design_candidates, args.n_reps)

    # Filter already-cached jobs in resume mode
    if args.resume:
        def is_cached(j):
            tag, rep = j[0]['tag'], j[1]
            out_json = os.path.join(args.work_dir, f"{tag}_rep{rep:03d}.json")
            if not os.path.exists(out_json):
                return False
            try:
                with open(out_json) as fh:
                    return 'ddG_switch_4corner' in json.load(fh)
            except Exception:
                return False
        before = len(jobs)
        jobs = [j for j in jobs if not is_cached(j)]
        print(f"  Resume: skipped {before - len(jobs):,} cached jobs, "
              f"{len(jobs):,} remaining")

    # ------------------------------------------------------------------
    # Run
    # ------------------------------------------------------------------
    all_results: dict[str, list[dict]] = defaultdict(list)
    t0 = time.time()

    print(f"\n  Running {len(jobs):,} jobs with {args.workers} workers...")
    with Pool(args.workers) as pool:
        for i, (tag, rep, result) in enumerate(
            pool.imap_unordered(score_one_replicate, jobs), 1
        ):
            all_results[tag].append(result)
            if i % 500 == 0:
                elapsed = time.time() - t0
                print(f"  {i:>7,}/{len(jobs):,}  elapsed={elapsed/3600:.1f}h  "
                      f"eta={elapsed/i*(len(jobs)-i)/3600:.1f}h")

    # ------------------------------------------------------------------
    # Aggregate replicates
    # ------------------------------------------------------------------
    # Build WT distribution from all WT replicates
    wt_ddgs = []
    for wt in wt_candidates:
        tag = wt['tag']
        for r in all_results.get(tag, []):
            if 'ddG_switch_4corner' in r:
                wt_ddgs.append(r['ddG_switch_4corner'])

    print(f"\n  Aggregating replicates...")
    summary_rows = []
    for candidate in wt_candidates + design_candidates:
        tag = candidate['tag']
        reps = all_results.get(tag, [])
        row = aggregate_replicates(
            tag, reps, wt_ddgs, args.outlier_threshold
        )
        # Carry over metadata from candidate record
        for k in ('run', 'hypothesis', 'subspace', 'sequence',
                  'n_mutations', 'mutations', 'mech_direction',
                  'mech_score_diff_HP_minus_ST', 'deviation_score_diff'):
            if k not in row:
                row[k] = candidate.get(k, '')
        summary_rows.append(row)

    write_csv(args.out_csv, summary_rows)

    elapsed = time.time() - t0
    print()
    print("=" * 78)
    print("  04_rosetta COMPLETE")
    print("=" * 78)
    print(f"  Time:     {elapsed/3600:.1f}h")
    print(f"  Output:   {args.out_csv} ({len(summary_rows)} rows)")

    # Quick stat summary
    ddgs = [r['ddG_median'] for r in summary_rows
            if 'ddG_median' in r and r.get('run') != 'TEMPLATE']
    if ddgs:
        import statistics
        print(f"\n  ddG_switch distribution ({len(ddgs)} designs):")
        print(f"    median: {statistics.median(ddgs):.2f} REU")
        print(f"    std:    {statistics.stdev(ddgs):.2f} REU")
        n_sig = sum(1 for r in summary_rows if r.get('significance') == 'significant')
        print(f"    Bonferroni-significant: {n_sig}")

    print(f"\n  Next step:")
    print(f"    python scripts/08_select_final.py --rosetta-csv {args.out_csv}")


if __name__ == "__main__":
    main()
