#!/usr/bin/env python3
"""
02_filter.py — plausibility filter and geometry-based mechanism scoring.

Processes all FASTAs from stage 01 (~448k sequences) and outputs a CSV of
plausible candidates with bidirectional phosphate-contact scores.

Plausibility filter (see src/phosphoswitch_design/filtering.py):
    - Length = 59 aa
    - Position 30 = Y (phospho-Tyr; kinase recognition)
    - Position 59 = V (C-terminal anchor)
    - <= MAX_MUTATIONS vs WT template (default 12)
    - Composition: hydrophobic, leucine, charged, aromatic percentages
    - No homopolymer run of 5+ identical residues
    Pass rate: ~25-30% of generated sequences.

Mechanism scoring (see src/phosphoswitch_design/mechanism.py):
    - Parse real phosphate centroids from phos_HP and phos_ST PDBs
    - For each candidate: score contacts on BOTH backbones
    - Record: HP score, ST score, donor differential, regional breakdown
    - deviation_score_diff = (HP−ST)_design − (HP−ST)_WT
    - H1/H2 matches_hypothesis flag

Usage
-----
    python scripts/02_filter.py --help
    python scripts/02_filter.py                              # defaults
    python scripts/02_filter.py --relaxed                    # bypass composition filters
    python scripts/02_filter.py --max-mutations 8 --max-kr 12
"""

import argparse
import csv
import glob
import os
import sys
import time
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from phosphoswitch_design.mechanism import (
    parse_pdb,
    mechanism_score,
    classify_h1_h2,
)
from phosphoswitch_design.filtering import (
    passes_plausibility,
    TEMPLATE,
)
from phosphoswitch_design.io_utils import (
    parse_fasta,
    parse_mpnn_header,
    get_subspace_from_path,
    get_state_temp_from_path,
    write_csv,
)


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--input-dir", default="outputs/01_generated_sequences",
                    help="Output dir from stage 01")
    ap.add_argument("--out-csv", default="outputs/02_plausible_with_mechanism.csv",
                    help="Output CSV path")
    ap.add_argument("--hairpin-pdb",
                    default="output/phase1/stateA_phospho_pulled.pdb",
                    help="Hairpin backbone PDB with phosphate (phos_HP)")
    ap.add_argument("--straight-pdb",
                    default="output/phase1/stateA_phospho.pdb",
                    help="Straight backbone PDB with phosphate (phos_ST)")
    ap.add_argument("--max-per-track", type=int, default=None,
                    help="Limit sequences per track (testing only)")
    ap.add_argument("--max-mutations", type=int, default=12)
    ap.add_argument("--max-charged", type=int, default=26,
                    help="Max total charged residues (WT=22)")
    ap.add_argument("--max-kr", type=int, default=14,
                    help="Max K+R count (WT=10)")
    ap.add_argument("--max-de", type=int, default=14,
                    help="Max D+E count (WT=12)")
    ap.add_argument("--relaxed", action="store_true",
                    help="Bypass biophysical composition filters. Still "
                         "enforces: length=59, Y30, V59, max_mutations, "
                         "no homopolymer>=5.")
    return ap


def main() -> None:
    ap = build_parser()
    args = ap.parse_args()

    os.makedirs(os.path.dirname(os.path.abspath(args.out_csv)), exist_ok=True)

    print("=" * 78)
    print("  02_filter — plausibility + mechanism scoring")
    print("=" * 78)
    if args.relaxed:
        print("\n  RELAXED MODE: composition filters bypassed.")
        print("  Enforced: length=59, Y30, V59, max_mutations, no homopolymer>=5")

    # ------------------------------------------------------------------
    # Parse reference backbone PDBs
    # ------------------------------------------------------------------
    for pdb_path, label in [(args.hairpin_pdb, "hairpin"), (args.straight_pdb, "straight")]:
        if not os.path.exists(pdb_path):
            sys.exit(f"ERROR: {label} PDB not found: {pdb_path}\n"
                     f"Set --hairpin-pdb and --straight-pdb to your phase1 PDB paths.")

    print(f"\n  Parsing backbone PDBs...")
    parsed_HP = parse_pdb(args.hairpin_pdb)
    parsed_ST = parse_pdb(args.straight_pdb)
    print(f"    Hairpin phos centroid:  {parsed_HP[0]}")
    print(f"    Straight phos centroid: {parsed_ST[0]}")

    # WT baseline mechanism
    wt_mech = mechanism_score(TEMPLATE, parsed_HP, parsed_ST)
    print(f"\n  WT baseline:")
    print(f"    HP donors={wt_mech['mech_HP_donors']}, "
          f"hbonds={wt_mech['mech_HP_hbonds']}, "
          f"repellers={wt_mech['mech_HP_repellers']}")
    print(f"    ST donors={wt_mech['mech_ST_donors']}, "
          f"hbonds={wt_mech['mech_ST_hbonds']}, "
          f"repellers={wt_mech['mech_ST_repellers']}")

    # ------------------------------------------------------------------
    # Discover tracks
    # ------------------------------------------------------------------
    if not os.path.isdir(args.input_dir):
        sys.exit(f"ERROR: input dir not found: {args.input_dir}\n"
                 f"Run 01_generate.py first.")

    tracks = []
    for d in sorted(os.listdir(args.input_dir)):
        full = os.path.join(args.input_dir, d)
        phase2_dir = os.path.join(full, "phase2")
        if os.path.isdir(phase2_dir):
            parts = d.split('_', 1)
            track_id = parts[0]
            label = parts[1] if len(parts) > 1 else d
            hyp = "H1" if track_id in ("1A", "2A") else ("H2" if track_id in ("1B", "2B") else "?")
            tracks.append((track_id, label, hyp, phase2_dir))

    print(f"\n  Found {len(tracks)} tracks: {[t[0] for t in tracks]}")

    # ------------------------------------------------------------------
    # Process all tracks
    # ------------------------------------------------------------------
    seen_seqs: dict[str, dict] = {}
    rejection_counts: Counter = Counter()
    total_scanned = 0
    t_start = time.time()

    for track_id, track_label, hypothesis, phase2_dir in tracks:
        fa_files = sorted(glob.glob(f"{phase2_dir}/**/*.fa", recursive=True))
        print(f"\n  Track {track_id} ({hypothesis}) — {len(fa_files)} FASTAs")

        n_total = n_passed = track_seqs = 0

        for fa in fa_files:
            subspace = get_subspace_from_path(fa)
            state, temp = get_state_temp_from_path(fa)

            for hdr, seq in parse_fasta(fa):
                n_total += 1
                if args.max_per_track and track_seqs >= args.max_per_track:
                    break

                ok, reason = passes_plausibility(
                    seq,
                    max_mutations=args.max_mutations,
                    max_total_charged=args.max_charged,
                    max_total_kr=args.max_kr,
                    max_total_de=args.max_de,
                    relaxed=args.relaxed,
                )
                if not ok:
                    rejection_counts[reason.split('=')[0]] += 1
                    continue

                n_passed += 1
                track_seqs += 1

                meta = parse_mpnn_header(hdr)
                muts = [(j+1, t, s) for j, (t, s) in enumerate(zip(TEMPLATE, seq)) if t != s]

                n_K = seq.count('K'); n_R = seq.count('R')
                n_D = seq.count('D'); n_E = seq.count('E')

                mech = mechanism_score(seq, parsed_HP, parsed_ST)

                hp_changed = mech['mech_HP_donors'] - wt_mech['mech_HP_donors']
                st_changed = mech['mech_ST_donors'] - wt_mech['mech_ST_donors']

                record = {
                    'run': track_id,
                    'hypothesis': hypothesis,
                    'subspace': subspace,
                    'mpnn_state': state,
                    'temperature': temp,
                    'sequence': seq,
                    'n_mutations': len(muts),
                    'mutations': ' '.join([f"{t}{p}{s}" for p, t, s in muts]),
                    'n_K': n_K, 'n_R': n_R, 'n_D': n_D, 'n_E': n_E,
                    'n_KR': n_K + n_R, 'n_DE': n_D + n_E,
                    'n_charged': n_K + n_R + n_D + n_E,
                    'mpnn_overall_conf': float(meta.get('overall_confidence', 0)),
                    'mpnn_ligand_conf': float(meta.get('ligand_confidence', 0)),
                    **mech,
                    'deviation_HP_donors': hp_changed,
                    'deviation_ST_donors': st_changed,
                    'deviation_HP_score': round(mech['mech_HP_score'] - wt_mech['mech_HP_score'], 2),
                    'deviation_ST_score': round(mech['mech_ST_score'] - wt_mech['mech_ST_score'], 2),
                    'deviation_score_diff': round(
                        mech['mech_score_diff_HP_minus_ST'] - wt_mech['mech_score_diff_HP_minus_ST'], 2),
                    'deviation_donor_diff': mech['mech_donor_diff_HP_minus_ST'] - wt_mech['mech_donor_diff_HP_minus_ST'],
                    'hp_donors_changed': hp_changed,
                    'st_donors_changed': st_changed,
                    'matches_hypothesis': classify_h1_h2(hypothesis, hp_changed, st_changed),
                }

                if seq not in seen_seqs or record['mpnn_overall_conf'] > seen_seqs[seq]['mpnn_overall_conf']:
                    seen_seqs[seq] = record

            if args.max_per_track and track_seqs >= args.max_per_track:
                break

        total_scanned += n_total
        print(f"    Total: {n_total:,}, Passed: {n_passed:,}")

    # WT reference row
    wt_record = {
        'run': 'TEMPLATE', 'hypothesis': 'WT', 'subspace': 'natural',
        'mpnn_state': 'WT', 'temperature': 'WT',
        'sequence': TEMPLATE, 'n_mutations': 0, 'mutations': '',
        'n_K': TEMPLATE.count('K'), 'n_R': TEMPLATE.count('R'),
        'n_D': TEMPLATE.count('D'), 'n_E': TEMPLATE.count('E'),
        'n_KR': TEMPLATE.count('K') + TEMPLATE.count('R'),
        'n_DE': TEMPLATE.count('D') + TEMPLATE.count('E'),
        'n_charged': sum(TEMPLATE.count(c) for c in 'KRDE'),
        'mpnn_overall_conf': 0.0, 'mpnn_ligand_conf': 0.0,
        **wt_mech,
        'deviation_HP_donors': 0, 'deviation_ST_donors': 0,
        'deviation_HP_score': 0, 'deviation_ST_score': 0,
        'deviation_score_diff': 0, 'deviation_donor_diff': 0,
        'hp_donors_changed': 0, 'st_donors_changed': 0,
        'matches_hypothesis': False,
    }

    all_kept = list(seen_seqs.values())
    all_kept.append(wt_record)

    # ------------------------------------------------------------------
    # Report
    # ------------------------------------------------------------------
    elapsed = time.time() - t_start
    print()
    print("=" * 78)
    print("  RESULTS")
    print("=" * 78)
    print(f"  Total scanned:    {total_scanned:,}")
    print(f"  Passed (unique):  {len(all_kept)-1:,}")
    print(f"  Time:             {elapsed/60:.1f} min")
    print(f"\n  Top rejection reasons:")
    for reason, count in rejection_counts.most_common(10):
        pct = count / total_scanned * 100 if total_scanned else 0
        print(f"    {reason:<20} {count:>10,}  ({pct:.1f}%)")

    direction_counts = Counter(r['mech_direction'] for r in all_kept if r['run'] != 'TEMPLATE')
    print(f"\n  Mechanism direction distribution:")
    for d, n in direction_counts.most_common():
        print(f"    {d:<25} {n:>10,}")

    write_csv(args.out_csv, all_kept)
    print(f"\n  Output: {args.out_csv} ({len(all_kept):,} rows)")
    print(f"\n  Next step:")
    print(f"    python scripts/03_select.py --input-csv {args.out_csv}")


if __name__ == "__main__":
    main()
