#!/usr/bin/env bash
# ============================================================================
# run_full_pipeline.sh — end-to-end phosphoswitch design pipeline
# ============================================================================
#
# Runs stages 01 → 02 → 03 → 04 → 08 in sequence for the LMNA Y45 target.
# Each stage checks that its required input exists before proceeding.
#
# PREREQUISITES
#   - conda activate phosphoswitch  (or equivalent env with all deps)
#   - LigandMPNN installed and path set in config/lmna_y45.yaml
#   - PyRosetta licensed and importable
#   - phase1 PDB files in output/phase1/ (stateA_phospho.pdb, etc.)
#
# USAGE
#   bash workflows/run_full_pipeline.sh                     # full run
#   bash workflows/run_full_pipeline.sh --resume            # resume from checkpoint
#   bash workflows/run_full_pipeline.sh --dry-run           # print plan, don't run
#   bash workflows/run_full_pipeline.sh --tracks 1A 1B      # subset of tracks
#   WORKERS=4 bash workflows/run_full_pipeline.sh           # custom parallelism
#   N_REPS=5   bash workflows/run_full_pipeline.sh          # quick validation test
#
# ESTIMATED COMPUTE (full run, default settings)
#   Stage 01 (LigandMPNN):  ~40h GPU    (4 tracks × 7 subspaces × 6 temps × 2 states)
#   Stage 02 (filter):      ~30 min CPU (single process, ~1.68M sequences)
#   Stage 03 (select):      <1 min CPU
#   Stage 04 (Rosetta):     ~25h CPU    (10k candidates × N=20 reps, 16 workers)
#   Stage 08 (consensus):   <1 min CPU
# ============================================================================

set -euo pipefail

# ── Configuration ────────────────────────────────────────────────────────────
CONFIG="${CONFIG:-config/lmna_y45.yaml}"
OUT_BASE="${OUT_BASE:-outputs}"
WORKERS="${WORKERS:-8}"
N_REPS="${N_REPS:-20}"
N_REPS_WT="${N_REPS_WT:-50}"
TRACKS="${TRACKS:-}"           # empty = all 4 tracks
RESUME="${RESUME:-}"           # non-empty = --resume flag
DRY_RUN="${DRY_RUN:-}"         # non-empty = --dry-run flag

# Parse command-line arguments
for arg in "$@"; do
    case "$arg" in
        --resume)   RESUME=1 ;;
        --dry-run)  DRY_RUN=1 ;;
        --tracks)   shift ;;   # handled below
        *)          ;;
    esac
done

# Build optional flags
RESUME_FLAG=""
if [ -n "$RESUME" ]; then RESUME_FLAG="--resume"; fi

DRY_FLAG=""
if [ -n "$DRY_RUN" ]; then DRY_FLAG="--dry-run"; fi

TRACKS_FLAG=""
if [ -n "$TRACKS" ]; then TRACKS_FLAG="--tracks $TRACKS"; fi

# Script directory (so we can run from any cwd)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

cd "$REPO_ROOT"

echo "============================================================"
echo "  phosphoswitch-design — full pipeline"
echo "============================================================"
echo "  Config:    $CONFIG"
echo "  Output:    $OUT_BASE"
echo "  Workers:   $WORKERS"
echo "  N reps:    $N_REPS (WT: $N_REPS_WT)"
echo "  Resume:    ${RESUME:-no}"
echo "  Dry run:   ${DRY_RUN:-no}"
echo ""

# ── Stage 01: LigandMPNN generation ─────────────────────────────────────────
echo "────────────────────────────────────────────────"
echo "  STAGE 01 — LigandMPNN generation"
echo "────────────────────────────────────────────────"
python scripts/01_generate.py \
    --config "$CONFIG" \
    --out-base "$OUT_BASE/01_generated_sequences" \
    $TRACKS_FLAG \
    $RESUME_FLAG \
    $DRY_FLAG

if [ -n "$DRY_RUN" ]; then
    echo "(Dry run: stopping after stage 01 plan)"
    exit 0
fi

# ── Stage 02: plausibility filter + mechanism scoring ────────────────────────
echo ""
echo "────────────────────────────────────────────────"
echo "  STAGE 02 — filter + mechanism score"
echo "────────────────────────────────────────────────"

# Use backbone PDB paths from default phase1 location.
# Override with --hairpin-pdb / --straight-pdb if your paths differ.
python scripts/02_filter.py \
    --input-dir  "$OUT_BASE/01_generated_sequences" \
    --out-csv    "$OUT_BASE/02_plausible_with_mechanism.csv" \
    --relaxed

# ── Stage 03: diversity selection ────────────────────────────────────────────
echo ""
echo "────────────────────────────────────────────────"
echo "  STAGE 03 — diverse top selection"
echo "────────────────────────────────────────────────"
python scripts/03_select.py \
    --input-csv "$OUT_BASE/02_plausible_with_mechanism.csv" \
    --out-csv   "$OUT_BASE/03_top_diverse_candidates.csv"

# ── Stage 04: deep Rosetta 4-state validation ────────────────────────────────
echo ""
echo "────────────────────────────────────────────────"
echo "  STAGE 04 — Rosetta 4-state validation"
echo "────────────────────────────────────────────────"
python scripts/04_rosetta.py \
    --input-csv  "$OUT_BASE/03_top_diverse_candidates.csv" \
    --out-csv    "$OUT_BASE/04_deep_validated.csv" \
    --work-dir   "$OUT_BASE/04_rosetta_cache" \
    --n-reps     "$N_REPS" \
    --n-reps-wt  "$N_REPS_WT" \
    --workers    "$WORKERS" \
    $RESUME_FLAG

# ── Stage 08: consensus ranking + final selection ────────────────────────────
echo ""
echo "────────────────────────────────────────────────"
echo "  STAGE 08 — consensus + final selection"
echo "────────────────────────────────────────────────"
python scripts/08_select_final.py \
    --rosetta-csv   "$OUT_BASE/04_deep_validated.csv" \
    --folding-csv   "$OUT_BASE/05_folding_solubility.csv" \
    --colabfold-csv "$OUT_BASE/06_colabfold_results.csv" \
    --af3-csv       "$OUT_BASE/07_af3_results.csv" \
    --out-dir       "$OUT_BASE/08_FINAL_WETLAB_candidates"

echo ""
echo "============================================================"
echo "  PIPELINE COMPLETE"
echo "============================================================"
echo ""
echo "  Final candidates: $OUT_BASE/08_FINAL_WETLAB_candidates/"
echo "    final_ranking.csv"
echo "    protein_sequences.fa"
echo "    codon_optimized_dna.fa"
echo "    selection_rationale.txt"
echo ""
echo "  Wet-lab next steps:"
echo "    1. Order WT + top candidates as 15N-labeled peptides"
echo "    2. Phosphorylate with Src kinase ± ATP"
echo "    3. Record NMR HSQC in ±phospho conditions"
echo "    4. Compare hairpin/straight conformational populations"
