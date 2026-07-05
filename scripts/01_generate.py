#!/usr/bin/env python3
"""
01_generate.py — LigandMPNN sequence generation for phosphoswitch design.

Calls LigandMPNN's run.py for all 4 hypothesis tracks × N subspaces ×
N temperatures × 2 backbone states (phos + apo).

Tracks
------
    1A  H1, straight+phos vs hairpin-apo    (phos stabilises STRAIGHT)
    1B  H2, hairpin+phos  vs straight-apo   (phos stabilises HAIRPIN)
    2A  H1, straight+phos vs B-hairpin-apo
    2B  H2, B-hairpin+phos vs straight-apo

One GPU process runs at a time.  No PyRosetta is used in this stage.

Usage
-----
    python scripts/01_generate.py --help
    python scripts/01_generate.py --config config/lmna_y45.yaml
    python scripts/01_generate.py --tracks 1A 1B --num-seqs 500 --dry-run
    python scripts/01_generate.py --resume   # skip already-completed runs
"""

import argparse
import os
import sys
import time
import yaml
from pathlib import Path

# Make the src package importable when running as a standalone script
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from phosphoswitch_design.sequence_gen import (
    TRACKS,
    run_ligandmpnn,
    prepare_track_pdbs,
    output_already_done,
)


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--config", default="config/lmna_y45.yaml",
                    help="Path to YAML config (default: config/lmna_y45.yaml)")
    ap.add_argument("--out-base", default="outputs/01_generated_sequences",
                    help="Output base directory")
    ap.add_argument("--num-seqs", type=int, default=None,
                    help="Sequences per (track, subspace, temp, state). "
                         "Overrides config value.")
    ap.add_argument("--temperatures", type=float, nargs="+", default=None,
                    help="Sampling temperatures (overrides config)")
    ap.add_argument("--tracks", nargs="+", default=None,
                    choices=["1A", "1B", "2A", "2B"],
                    help="Subset of tracks to run (default: all)")
    ap.add_argument("--subspaces", nargs="+", default=None,
                    help="Subset of subspace names to run (default: all)")
    ap.add_argument("--seed", type=int, default=42,
                    help="LigandMPNN random seed")
    ap.add_argument("--resume", action="store_true",
                    help="Skip (track, subspace, temp, state) combos that "
                         "already produced FASTAs")
    ap.add_argument("--dry-run", action="store_true",
                    help="Print what would run without executing")
    return ap


def main() -> None:
    ap = build_parser()
    args = ap.parse_args()

    # ------------------------------------------------------------------
    # Load config
    # ------------------------------------------------------------------
    if not os.path.exists(args.config):
        sys.exit(f"ERROR: config not found: {args.config}")

    with open(args.config) as fh:
        cfg = yaml.safe_load(fh)

    ligandmpnn_dir  = cfg["ligand_mpnn_dir"]
    ligandmpnn_ckpt = cfg["ligand_mpnn_checkpoint"]
    num_seqs        = args.num_seqs or cfg.get("mpnn_num_sequences", 2000)
    temperatures    = args.temperatures or cfg.get("mpnn_temperatures", [0.1, 0.2, 0.3, 0.5])
    omit_aa         = cfg.get("omit_aa_design") or []
    template_seq    = cfg.get("template_seq",
                              "ASSTPLSPTRITRLQEKEDLQELNDRLAVYIDRVRSLETENAGLRLRITESEEVVSREV")
    n_residues      = len(template_seq)

    all_subspaces = cfg.get("design_subspaces", [])
    if not all_subspaces:
        sys.exit("ERROR: no design_subspaces in config")

    # ------------------------------------------------------------------
    # Filter subspaces / tracks
    # ------------------------------------------------------------------
    if args.subspaces:
        wanted = set(args.subspaces)
        subspaces = [s for s in all_subspaces if s["name"] in wanted]
        if not subspaces:
            avail = [s["name"] for s in all_subspaces]
            sys.exit(f"ERROR: none of {args.subspaces} match available subspaces: {avail}")
    else:
        subspaces = all_subspaces

    if args.tracks:
        wanted_tracks = set(args.tracks)
        tracks = [t for t in TRACKS if t["id"] in wanted_tracks]
    else:
        tracks = TRACKS

    # ------------------------------------------------------------------
    # Sanity checks
    # ------------------------------------------------------------------
    if not os.path.isfile(os.path.join(ligandmpnn_dir, "run.py")):
        sys.exit(f"ERROR: LigandMPNN run.py not found in {ligandmpnn_dir}")
    if not os.path.isfile(ligandmpnn_ckpt):
        sys.exit(f"ERROR: checkpoint not found: {ligandmpnn_ckpt}")

    # Verify PDBs; skip tracks with missing files
    valid_tracks = []
    for t in tracks:
        phos_ok = os.path.isfile(t["phos_pdb"])
        apo_ok  = os.path.isfile(t["apo_pdb"])
        if phos_ok and apo_ok:
            valid_tracks.append(t)
        else:
            print(f"  [SKIP] Track {t['id']} — missing PDB:")
            if not phos_ok:
                print(f"    phos: {t['phos_pdb']}")
            if not apo_ok:
                print(f"    apo:  {t['apo_pdb']}")

    if not valid_tracks:
        sys.exit("ERROR: no tracks have both PDBs present")

    # ------------------------------------------------------------------
    # Print plan
    # ------------------------------------------------------------------
    n_combos = len(valid_tracks) * len(subspaces) * len(temperatures) * 2
    total_seqs = n_combos * num_seqs

    print("=" * 78)
    print("  01_generate — v2 LigandMPNN sequence generation")
    print("=" * 78)
    print(f"  Tracks ({len(valid_tracks)}):       {[t['id'] for t in valid_tracks]}")
    print(f"  Subspaces ({len(subspaces)}):     {[s['name'] for s in subspaces]}")
    print(f"  Temperatures:        {temperatures}")
    print(f"  Seqs per combo:      {num_seqs}")
    print(f"  Total MPNN runs:     {n_combos}")
    print(f"  Total seqs:          ~{total_seqs:,}")
    print(f"  Output base:         {args.out_base}")
    print(f"  Resume:              {args.resume}")
    print(f"  Omit AA:             {omit_aa or '(none)'}")
    print(f"  LigandMPNN dir:      {ligandmpnn_dir}")
    print(f"  Est GPU time:        ~{n_combos * num_seqs * 0.2 / 3600:.1f}h")

    if args.dry_run:
        print("\n  --dry-run: exiting without executing.")
        return

    # ------------------------------------------------------------------
    # Run
    # ------------------------------------------------------------------
    os.makedirs(args.out_base, exist_ok=True)
    overall_start = time.time()
    n_done = n_skipped = n_failed = 0

    for ti, track in enumerate(valid_tracks, 1):
        track_dir = Path(args.out_base) / f"{track['id']}_{track['label']}"
        phase2_root = track_dir / "phase2"
        phase2_root.mkdir(parents=True, exist_ok=True)

        phos_with_ligand, apo_clean = prepare_track_pdbs(track, phase2_root)

        print()
        print("=" * 78)
        print(f"  TRACK {ti}/{len(valid_tracks)}: {track['id']} ({track['hypothesis']}) "
              f"— {track['label']}")
        print(f"    phos → {phos_with_ligand}")
        print(f"    apo  → {apo_clean}")
        print("=" * 78)

        for si, subspace in enumerate(subspaces, 1):
            name      = subspace["name"]
            positions = subspace["positions"]
            fixed     = [r for r in range(1, n_residues + 1) if r not in positions]

            sub_dir = phase2_root / name
            sub_dir.mkdir(exist_ok=True)

            print(f"\n  [{si}/{len(subspaces)}] {name}  ({len(positions)} designable)")

            for state_label, mpnn_pdb, has_ligand in [
                ("A", str(phos_with_ligand), True),
                ("B", str(apo_clean),        False),
            ]:
                state_dir = sub_dir / f"lmpnn_out_{state_label}"
                state_dir.mkdir(exist_ok=True)

                for temp in temperatures:
                    temp_dir = state_dir / f"T{temp}"

                    if args.resume and output_already_done(str(temp_dir)):
                        n_skipped += 1
                        continue

                    t0 = time.time()
                    success, msg = run_ligandmpnn(
                        pdb_path=mpnn_pdb,
                        out_dir=str(temp_dir),
                        has_ligand=has_ligand,
                        design_positions=positions,
                        fixed_residues=fixed,
                        num_seqs=num_seqs,
                        temperature=temp,
                        ligandmpnn_dir=ligandmpnn_dir,
                        ligandmpnn_ckpt=ligandmpnn_ckpt,
                        seed=args.seed,
                        omit_aa=omit_aa,
                    )
                    dt = time.time() - t0

                    if success:
                        n_done += 1
                        n_fa = len(list((temp_dir / "seqs").glob("*.fa"))) \
                               if (temp_dir / "seqs").is_dir() else 0
                        print(f"      state={state_label} T={temp}  "
                              f"ok ({dt:.0f}s, {n_fa} fastas)")
                    else:
                        n_failed += 1
                        print(f"      state={state_label} T={temp}  "
                              f"FAILED: {msg[:200]}")

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    elapsed = time.time() - overall_start
    print()
    print("=" * 78)
    print("  GENERATION COMPLETE")
    print("=" * 78)
    print(f"  Total time:   {elapsed/3600:.1f}h")
    print(f"  Successful:   {n_done}")
    print(f"  Failed:       {n_failed}")
    print(f"  Skipped:      {n_skipped}")
    print(f"\n  Next step:")
    print(f"    python scripts/02_filter.py --input-dir {args.out_base}")


if __name__ == "__main__":
    main()
