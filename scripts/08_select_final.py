#!/usr/bin/env python3
"""
08_select_final.py — consensus ranking and final wet-lab candidate selection.

Combines all signals from previous pipeline stages into a composite consensus
score, then selects 8-12 candidates with diversity constraints.

Input CSVs (all optional except --rosetta-csv)
-----------------------------------------------
    04_deep_validated.csv       Rosetta z-score, ddG statistics
    05_folding_solubility.csv   ESMFold/OmegaFold stability + solubility
    06_colabfold_results.csv    ColabFold 1.5.5 AF2 bistability
    07_af3_results.csv          AF3 3.0.3 with PTR ligand, pLDDT, RMSD

Consensus scoring weights (per candidate)
-----------------------------------------
    Rosetta z-score     ×2.0 per |z| unit   primary thermodynamic signal
    Mechanism dev       ×0.15 per unit       geometry contact differential
    Direction agree     1-3 pts              Rosetta+mech+AF direction match
    Folding stability   1.5 pts              apo energies < 0 REU
    Solubility          0.5-1.0 pts          predicted solubility score
    AF2 bistability     0.5-1.0 pts          apo↔phos RMSD conformational signal
    AF3 confidence      1.0 pts              apo + phos pLDDT > 70
    AF3 change          1.5 pts              AF3 predicts conformational change
    Minimal mutations   0.5-1.0 pts          ≤7 muts / ≤5 muts bonus
    No outlier          0.5 pts              Rosetta range ≤ 15 REU

Penalties:
    H-cluster           -5.0                 Always (PyRosetta artifact)
    High aggregation    -2.0                 Strict mode (ΔAGG vs WT > 0.10)
    PTM crosstalk       -1.0 each            Strict mode

Hard gate: passes_physical_filters must be True.

Outputs
-------
    final_ranking.csv           ranked table with all scores
    protein_sequences.fa        FASTA of final set + WT control
    codon_optimized_dna.fa      E. coli optimised DNA (for synthesis order)
    selection_rationale.txt     human-readable decision log

Usage
-----
    python scripts/08_select_final.py --help
    python scripts/08_select_final.py --rosetta-csv outputs/04_deep_validated.csv
    python scripts/08_select_final.py --n-final 8 --relaxed
"""

import argparse
import csv
import os
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from phosphoswitch_design.consensus import (
    compute_consensus,
    select_final_candidates,
    reverse_translate,
    get_mutations,
    TEMPLATE,
)
from phosphoswitch_design.filtering import get_mutation_signature
from phosphoswitch_design.io_utils import load_csv_indexed, write_csv


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--rosetta-csv",    default="outputs/04_deep_validated.csv")
    ap.add_argument("--folding-csv",    default="outputs/05_folding_solubility.csv")
    ap.add_argument("--colabfold-csv",  default="outputs/06_colabfold_results.csv")
    ap.add_argument("--af3-csv",        default="outputs/07_af3_results.csv")
    ap.add_argument("--out-dir",        default="outputs/08_FINAL_WETLAB_candidates")
    ap.add_argument("--n-final",        type=int, default=12,
                    help="Number of final candidates to select (default 12)")
    ap.add_argument("--min-directions", type=int, default=2,
                    help="Minimum distinct mechanism directions in selection")
    ap.add_argument("--min-signatures", type=int, default=4,
                    help="Minimum distinct mutation signatures in selection")
    ap.add_argument("--relaxed",        action="store_true",
                    help="Skip aggregation and PTM crosstalk penalties")
    return ap


def main() -> None:
    ap = build_parser()
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    print("=" * 78)
    print("  08_select_final — consensus ranking + wet-lab selection")
    print("=" * 78)
    if args.relaxed:
        print("  RELAXED MODE: aggregation + PTM penalties skipped")

    # ------------------------------------------------------------------
    # Load all CSVs
    # ------------------------------------------------------------------
    print(f"\n  Loading data sources...")
    rosetta_data   = load_csv_indexed(args.rosetta_csv)
    folding_data   = load_csv_indexed(args.folding_csv)
    colabfold_data = load_csv_indexed(args.colabfold_csv)
    af3_data       = load_csv_indexed(args.af3_csv)

    print(f"    Rosetta:    {len(rosetta_data)}")
    print(f"    Folding:    {len(folding_data)}")
    print(f"    ColabFold:  {len(colabfold_data)}")
    print(f"    AF3:        {len(af3_data)}")

    if not rosetta_data:
        sys.exit(f"ERROR: no data in {args.rosetta_csv}\nRun 04_rosetta.py first.")

    # ------------------------------------------------------------------
    # Merge all signals
    # ------------------------------------------------------------------
    all_tags = set(rosetta_data.keys()) | set(folding_data.keys())
    merged = []
    for tag in all_tags:
        rec: dict = {}
        rec.update(folding_data.get(tag, {}))
        rec.update(rosetta_data.get(tag, {}))
        # ColabFold and AF3: only fill in missing keys
        for k, v in colabfold_data.get(tag, {}).items():
            if v and v not in ('None', ''):
                rec.setdefault(k, v)
        for k, v in af3_data.get(tag, {}).items():
            if v and v not in ('None', ''):
                rec.setdefault(k, v)
        rec['tag'] = tag
        merged.append(rec)

    print(f"\n  Total merged: {len(merged):,} candidates")

    # ------------------------------------------------------------------
    # Compute consensus scores
    # ------------------------------------------------------------------
    print(f"  Computing consensus scores...")
    for rec in merged:
        rec['consensus_score'], rec['consensus_contribs'] = compute_consensus(
            rec, relaxed=args.relaxed
        )
        rec['mutation_signature'] = get_mutation_signature(rec.get('sequence', ''))

    # ------------------------------------------------------------------
    # Apply hard gate and select
    # ------------------------------------------------------------------
    wt_rec   = next((r for r in merged if r.get('tag') == 'WT_TEMPLATE'), None)
    designs  = [r for r in merged if r.get('tag') != 'WT_TEMPLATE']

    # Hard gate: passes_physical_filters
    designs = [r for r in designs
               if str(r.get('passes_physical_filters', '')).lower() == 'true']
    print(f"  After physical filter: {len(designs):,}")

    selected = select_final_candidates(
        designs,
        n_final=args.n_final,
        min_directions=args.min_directions,
        min_signatures=args.min_signatures,
    )

    dirs_seen = {r.get('mech_direction', '?') for r in selected}
    sigs_seen = {r.get('mutation_signature', '?') for r in selected}
    print(f"\n  Selected: {len(selected)}")
    print(f"  Directions: {dirs_seen}")
    print(f"  Signatures: {len(sigs_seen)}")

    # ------------------------------------------------------------------
    # Write outputs
    # ------------------------------------------------------------------

    # 1. Final ranking CSV
    out_csv = os.path.join(args.out_dir, "final_ranking.csv")
    rank_rows = []
    for rank, rec in enumerate(selected, 1):
        seq = rec.get('sequence', '')
        rank_rows.append({
            'rank': rank,
            'tag': rec.get('tag'),
            'sequence': seq,
            'mutations': ','.join(get_mutations(seq)),
            'n_mutations': rec.get('n_mutations'),
            'mutation_signature': rec.get('mutation_signature'),
            'orig_run': rec.get('orig_run', rec.get('run', '')),
            'orig_subspace': rec.get('orig_subspace', rec.get('subspace', '')),
            'orig_hypothesis': rec.get('orig_hypothesis', rec.get('hypothesis', '')),
            'consensus_score': rec.get('consensus_score'),
            'rosetta_z_score': rec.get('z_score_vs_WT'),
            'rosetta_significance': rec.get('significance'),
            'rosetta_ddG_median': rec.get('ddG_median'),
            'mech_direction': rec.get('mech_direction'),
            'mech_deviation': rec.get('mech_deviation', rec.get('deviation_score_diff', '')),
            'solubility_score': rec.get('solubility_score'),
            'af2_apo_vs_phos_rmsd': rec.get('apo_vs_phos_rmsd'),
            'af3_apo_plddt': rec.get('af3_apo_plddt'),
            'af3_predicts_change': rec.get('af3_predicts_change'),
            'consensus_contribs': str(rec.get('consensus_contribs', {})),
        })
    write_csv(out_csv, rank_rows)
    print(f"\n  {out_csv}")

    # 2. Protein FASTA
    out_fa = os.path.join(args.out_dir, "protein_sequences.fa")
    with open(out_fa, 'w') as fh:
        fh.write(f">WT_LMNA_phospho_switch | template | length={len(TEMPLATE)}\n"
                 f"{TEMPLATE}\n\n")
        for rank, rec in enumerate(selected, 1):
            seq  = rec.get('sequence', '')
            tag  = rec.get('tag', f'cand_{rank}')
            muts = ','.join(get_mutations(seq))
            sc   = rec.get('consensus_score', 0)
            z    = rec.get('z_score_vs_WT', '?')
            fh.write(f">rank{rank:02d}_{tag} | muts={muts} | score={sc} | z={z}\n"
                     f"{seq}\n\n")
    print(f"  {out_fa}")

    # 3. Codon-optimised DNA
    out_dna = os.path.join(args.out_dir, "codon_optimized_dna.fa")
    with open(out_dna, 'w') as fh:
        fh.write(f">WT_LMNA_phospho_switch_DNA | E.coli optimized\n"
                 f"{reverse_translate(TEMPLATE)}\n\n")
        for rank, rec in enumerate(selected, 1):
            seq = rec.get('sequence', '')
            tag = rec.get('tag', f'cand_{rank}')
            fh.write(f">rank{rank:02d}_{tag}_DNA | E.coli optimized\n"
                     f"{reverse_translate(seq)}\n\n")
    print(f"  {out_dna}")

    # 4. Selection rationale
    out_txt = os.path.join(args.out_dir, "selection_rationale.txt")
    with open(out_txt, 'w') as fh:
        fh.write("=" * 78 + "\n  WET-LAB CANDIDATE SELECTION RATIONALE\n"
                 + "=" * 78 + "\n\n")
        fh.write("CONSENSUS SCORING WEIGHTS:\n"
                 "  Rosetta z-score:        ×2.0 per |z|\n"
                 "  Mechanism deviation:    ×0.15 per unit (capped 10)\n"
                 "  Direction agreement:    1-3 pts\n"
                 "  Folding stability:      1.5 pts\n"
                 "  Solubility:             0.5-1.0 pts\n"
                 "  AF2 bistability:        0.5-1.0 pts\n"
                 "  AF3 confidence:         1.0 pts\n"
                 "  AF3 conformational chg: 1.5 pts\n"
                 "  Minimal mutations:      0.5-1.0 pts\n"
                 "  No outlier:             0.5 pts\n\n"
                 "PENALTIES:\n"
                 "  H-cluster:     -5.0 (always — PyRosetta artifact)\n"
                 "  High agg:      -2.0 (strict mode)\n"
                 "  PTM crosstalk: -1.0 each (strict mode)\n\n")
        fh.write("=" * 78 + "\n  SELECTED CANDIDATES\n" + "=" * 78 + "\n\n")
        for rank, rec in enumerate(selected, 1):
            seq  = rec.get('sequence', '')
            muts = ','.join(get_mutations(seq))
            fh.write(f"\nRANK {rank}: {rec.get('tag')}\n" + "-" * 60 + "\n")
            fh.write(f"  Sequence:      {seq}\n")
            fh.write(f"  Mutations:     {muts}\n")
            fh.write(f"  N mutations:   {rec.get('n_mutations')}\n")
            fh.write(f"  Signature:     {rec.get('mutation_signature')}\n")
            fh.write(f"\n  CONSENSUS SCORE: {rec.get('consensus_score')}\n")
            contribs = rec.get('consensus_contribs', {})
            if isinstance(contribs, dict):
                for k, v in contribs.items():
                    fh.write(f"    {k:<25} {v}\n")
            fh.write(f"\n  Rosetta z:     {rec.get('z_score_vs_WT')}\n")
            fh.write(f"  ddG median:    {rec.get('ddG_median')} REU\n")
            fh.write(f"  Direction:     {rec.get('mech_direction')}\n")
            fh.write(f"  AF2 RMSD:      {rec.get('apo_vs_phos_rmsd')}\n")
    print(f"  {out_txt}")

    # ------------------------------------------------------------------
    # Console summary
    # ------------------------------------------------------------------
    print()
    print("=" * 78)
    print("  FINAL SELECTION")
    print("=" * 78)
    print(f"\n  {'Rank':>4} {'Tag':<35} {'Score':>6} {'Z':>6} {'Dir':<14} {'Muts':>4}")
    print(f"  {'-'*4} {'-'*35} {'-'*6} {'-'*6} {'-'*14} {'-'*4}")
    for rank, rec in enumerate(selected, 1):
        z = rec.get('z_score_vs_WT', '?')
        try:
            z = f"{float(z):.2f}"
        except (TypeError, ValueError):
            pass
        print(f"  {rank:>4} {rec.get('tag', '?'):<35} "
              f"{rec.get('consensus_score', 0):>6.2f} {z:>6} "
              f"{rec.get('mech_direction', '?'):<14} "
              f"{rec.get('n_mutations', '?'):>4}")

    print(f"\n  Outputs: {args.out_dir}/")
    print(f"\n  Wet-lab plan:")
    print(f"    1. Order WT + {len(selected)} candidates as 15N-labeled peptides")
    print(f"    2. NMR HSQC ± phosphorylation by Src kinase on each")
    print(f"    3. Compare hairpin/straight peak populations between ±phos")


if __name__ == "__main__":
    main()
