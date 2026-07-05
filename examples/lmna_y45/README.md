# LMNA Y45 Example — Multi-State Design

Complete worked example of the phosphoswitch design pipeline applied to LMNA Y45 (Src kinase site).

## Backbone structures

| File | Description |
|------|-------------|
| `backbones/stateA_phospho.pdb` | Straight helix backbone, phosphotyrosine at position 30 (Y45 in full LMNA). Used by H1 tracks as the phospho target state. |
| `backbones/stateA_phospho_pulled.pdb` | Hairpin backbone, phosphotyrosine at position 30. Used by H2 tracks as the phospho target state. RMSD vs stateA_phospho = 34.35 Å — this is the switch magnitude. |
| `backbones/stateB_aln.pdb` | Natural hairpin backbone (apo). Used by tracks 2A/2B. |
| `backbones/stateA_aln.pdb` | Natural straight backbone (apo). |

## Top designed candidates

`top_candidates/boltz_ranked_clean_switches.csv` — 36 candidates that passed all filters:
- 4-state Rosetta ddG_switch < −5 REU (N=10 replicates)
- ColabFold RMSD between states > 5 Å
- Boltz 2.2.1 ranked structure: apo vs phos RMSD > 3 Å
- No histidine cluster artifact (≤5 H per candidate)

Top candidate: `cand_24213_2A_bend_plus_face` — Boltz RMSD 22.4 Å between apo and phos states.

## Running the pipeline

```bash
# Stage 01: generate sequences
python scripts/01_generate.py --config examples/lmna_y45/config_motif_fixed_v2.yaml \
    --out-base outputs/01_generated --num-seqs 5000

# Stage 02: filter by mechanism
python scripts/02_filter.py --input-dir outputs/01_generated --out-csv outputs/02_filtered.csv

# Stage 03: select diverse top candidates
python scripts/03_select.py --input-csv outputs/02_filtered.csv --out-csv outputs/03_top.csv

# Stage 04: Rosetta 4-state validation
python scripts/04_rosetta.py --input-csv outputs/03_top.csv --out-csv outputs/04_rosetta.csv \
    --phos-st examples/lmna_y45/backbones/stateA_phospho.pdb \
    --phos-hp examples/lmna_y45/backbones/stateA_phospho_pulled.pdb
```
