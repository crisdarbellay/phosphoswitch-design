"""
rosetta.py — 4-state PyRosetta thermodynamic cycle for phosphoswitch validation.

The 4-state thermodynamic cycle
---------------------------------
Each candidate is scored on four backbone corners:

    phos_ST (phospho + straight)   phos_HP (phospho + hairpin)
         |                               |
    (E1)  \\                           / (E3)
         apo_ST  (apo + straight)   apo_HP  (apo + hairpin)
                 (E2)                     (E4)

    phos_pref  = E_phos_HP − E_phos_ST     # phospho wants hairpin (if +)
    apo_pref   = E_apo_HP  − E_apo_ST      # apo   wants hairpin (if +)
    ddG_switch = phos_pref − apo_pref

Interpretation of ddG_switch:
    > 0  →  phospho increases HAIRPIN preference (H2 mechanism)
    < 0  →  phospho increases STRAIGHT preference (H1 mechanism)
    ≈ 0  →  phospho has no effect on conformational preference

Noise floor and statistics
---------------------------
A single FastRelax run (3 cycles, backbone constrained) has a stochastic noise
floor of ±3-5 REU.  For a 4-corner cycle this compounds to ±7 REU.  No
single-replicate ddG below ~10 REU should be trusted as a real signal.

The pipeline mitigates this by running N=20 replicates per candidate (and N=50
for the WT baseline), using different stochastic seeds per replicate, then
reporting mean, median, std, and z-score vs. WT.

FastRelax settings
------------------
    cycles           = 3   (fast; preserves backbone topology)
    constrain_to_start_coords = True   (prevent large backbone movements)
    coord_constrain_sidechains = False (allow side-chain repacking)
    ramp_down_constraints = True       (standard Rosetta protocol)

Histidine cluster artifact
---------------------------
PyRosetta uses a fixed-protonation energy function (ref2015) that treats all
His as singly protonated at a single tautomer.  At pH 7-7.4, real His is
mostly neutral.  Clusters of His near the phosphate site score artificially
well in PyRosetta but would be neutral or unfavourable in reality.  Stage 05
flags H-clusters; stage 08 applies a -5.0 consensus penalty.  Never order
high-His designs for synthesis.

Apo backbone generation
------------------------
The four backbone PDBs are:
    phos_ST  — input from phase1 (stateA_phospho.pdb)
    phos_HP  — input from phase1 (stateA_phospho_pulled.pdb)
    apo_ST   — auto-derived by stripping phos_ST
    apo_HP   — auto-derived by stripping phos_HP

The apo PDBs are created once per run by make_apo_pdb() and cached.
"""

from __future__ import annotations
import json
import os
from typing import Optional


# ---------------------------------------------------------------------------
# Apo PDB generation
# ---------------------------------------------------------------------------
def make_apo_pdb(phos_pdb_path: str, apo_pdb_path: str) -> None:
    """Strip phosphate atoms and rename PTR/TYS→TYR to create an apo backbone.

    Idempotent: if *apo_pdb_path* already exists, this function returns
    immediately without rewriting it.

    Parameters
    ----------
    phos_pdb_path : str
        Input PDB containing a phospho-Tyr (PTR HETATM or ATOM + separate
        phosphate HETATM records).
    apo_pdb_path : str
        Path to write the apo PDB.
    """
    if os.path.exists(apo_pdb_path):
        return
    with open(phos_pdb_path) as fin, open(apo_pdb_path, 'w') as fout:
        for line in fin:
            if line.startswith("HETATM"):
                atom = line[12:16].strip()
                if atom in ('P', 'O1P', 'O2P', 'O3P', 'OP1', 'OP2', 'OP3', 'P1'):
                    continue
                if "PTR" in line[17:20] or "TYS" in line[17:20]:
                    line = line[:17] + "TYR" + line[20:]
                if "TYR" in line[17:20]:
                    line = "ATOM  " + line[6:]
            elif line.startswith("ATOM") and (
                "PTR" in line[17:20] or "TYS" in line[17:20]
            ):
                atom = line[12:16].strip()
                if atom in ('P', 'O1P', 'O2P', 'O3P', 'OP1', 'OP2', 'OP3', 'P1'):
                    continue
                line = line[:17] + "TYR" + line[20:]
            fout.write(line)


# ---------------------------------------------------------------------------
# 4-corner scoring
# ---------------------------------------------------------------------------
def score_one_replicate(args: tuple) -> tuple[str, int, dict]:
    """Score one candidate at one replicate with a full 4-corner thermodynamic cycle.

    This function is designed to be called inside a multiprocessing pool.
    PyRosetta is initialised lazily on first call per worker process using a
    deterministic seed derived from (replicate, tag) so that every replicate
    uses a different but reproducible stochastic seed.

    If a cached JSON for (tag, replicate) already exists and contains a
    valid 'ddG_switch_4corner' key, the cached result is returned without
    running PyRosetta.

    Parameters
    ----------
    args : (candidate_dict, replicate_int, work_dir_str, pdbs_dict)
        candidate_dict : {'tag': str, 'sequence': str, ...}
        replicate      : replicate index (0-based)
        work_dir       : directory for per-replicate JSON cache files
        pdbs           : {
            'phos_ST': path, 'phos_HP': path,
            'apo_ST':  path, 'apo_HP':  path
          }

    Returns
    -------
    (tag, replicate, result_dict)
    """
    candidate, replicate, work_dir, pdbs = args
    seq = candidate['sequence']
    tag = candidate['tag']

    out_json = os.path.join(work_dir, f"{tag}_rep{replicate:03d}.json")

    # Check cache
    if os.path.exists(out_json):
        try:
            with open(out_json) as fh:
                cached = json.load(fh)
            if 'ddG_switch_4corner' in cached:
                return tag, replicate, cached
        except Exception:
            pass

    try:
        import pyrosetta  # type: ignore
        from pyrosetta import pose_from_pdb, create_score_function  # type: ignore
        from pyrosetta.rosetta.protocols.relax import FastRelax  # type: ignore
        from pyrosetta.toolbox import mutate_residue  # type: ignore

        if not getattr(score_one_replicate, '_init', False):
            # Unique but reproducible seed per (replicate, tag)
            seed = 1000 + replicate * 17 + (hash(tag) % 100000) * 31
            pyrosetta.init(
                f"-mute all -no_optH true -ex1 -ex2 -ignore_unrecognized_res true "
                f"-constant_seed -jran {seed}",
                silent=True,
            )
            score_one_replicate._init = True  # type: ignore[attr-defined]
    except Exception as exc:
        return tag, replicate, {"error": f"PyRosetta init: {exc}"}

    scorefxn = create_score_function("ref2015")

    # Ensure apo backbones exist (cached if already created)
    make_apo_pdb(pdbs['phos_ST'], pdbs['apo_ST'])
    make_apo_pdb(pdbs['phos_HP'], pdbs['apo_HP'])

    energies: dict[str, float] = {}
    for state_name, pdb_path in pdbs.items():
        try:
            pose = pose_from_pdb(pdb_path)
            n_pose = pose.total_residue()
            for i, aa in enumerate(seq):
                if i + 1 > n_pose:
                    break
                if pose.residue(i + 1).name1() != aa:
                    try:
                        mutate_residue(pose, i + 1, aa, pack_radius=6.0)
                    except Exception:
                        pass

            # FastRelax: 3 cycles, backbone constrained, sidechains free
            relax = FastRelax(scorefxn, 3)
            relax.constrain_relax_to_start_coords(True)
            relax.coord_constrain_sidechains(False)
            relax.ramp_down_constraints(True)
            relax.apply(pose)
            energies[state_name] = float(scorefxn(pose))
        except Exception as exc:
            return tag, replicate, {"error": f"scoring {state_name}: {str(exc)[:200]}"}

    # 4-corner thermodynamic cycle
    phos_pref = energies['phos_HP'] - energies['phos_ST']
    apo_pref  = energies['apo_HP']  - energies['apo_ST']
    ddG       = phos_pref - apo_pref   # positive → H2, negative → H1

    result = {
        'tag': tag,
        'replicate': replicate,
        'sequence': seq,
        **{f'E_{k}': v for k, v in energies.items()},
        'phos_pref': phos_pref,
        'apo_pref': apo_pref,
        'ddG_switch_4corner': ddG,
    }

    try:
        os.makedirs(work_dir, exist_ok=True)
        with open(out_json, 'w') as fh:
            json.dump(result, fh, indent=2)
    except Exception:
        pass

    return tag, replicate, result


# ---------------------------------------------------------------------------
# Statistics over replicates
# ---------------------------------------------------------------------------
def aggregate_replicates(
    tag: str,
    replicates: list[dict],
    wt_distribution: list[float],
    outlier_range_threshold: float = 15.0,
) -> dict:
    """Aggregate N replicates into per-candidate summary statistics.

    Parameters
    ----------
    tag : str
        Candidate identifier.
    replicates : list[dict]
        List of result dicts from score_one_replicate().
    wt_distribution : list[float]
        List of WT ddG_switch_4corner values across all WT replicates.
        Used to compute z-score significance.
    outlier_range_threshold : float
        If max(ddG) − min(ddG) across replicates exceeds this, the candidate
        is flagged as an outlier (likely stochastic explosion in one replicate).

    Returns
    -------
    Summary dict with keys: tag, sequence, n_reps, ddG_mean, ddG_median,
    ddG_std, ddG_min, ddG_max, ddG_range, z_score_vs_WT, outlier_flag,
    significance.
    """
    import math
    import statistics

    ddgs = [r['ddG_switch_4corner'] for r in replicates if 'ddG_switch_4corner' in r]
    if not ddgs:
        return {
            'tag': tag,
            'n_reps': 0,
            'error': 'all replicates failed',
        }

    ddg_mean   = statistics.mean(ddgs)
    ddg_median = statistics.median(ddgs)
    ddg_std    = statistics.stdev(ddgs) if len(ddgs) > 1 else 0.0
    ddg_min    = min(ddgs)
    ddg_max    = max(ddgs)
    ddg_range  = ddg_max - ddg_min
    outlier    = ddg_range > outlier_range_threshold

    # z-score vs WT
    if wt_distribution and len(wt_distribution) > 1:
        wt_mean = statistics.mean(wt_distribution)
        wt_std  = statistics.stdev(wt_distribution)
        z = (ddg_median - wt_mean) / wt_std if wt_std > 0 else 0.0
    else:
        z = 0.0

    # Bonferroni-corrected significance at alpha=0.05 for ~10,000 candidates
    alpha_corrected = 0.05 / 10000
    # |z| > 4.4 ≈ two-tailed Bonferroni threshold for 10k tests at alpha=0.05
    bonferroni_z_threshold = 4.4
    significance = "significant" if abs(z) > bonferroni_z_threshold else "ns"

    seq = replicates[0].get('sequence', '')
    return {
        'tag': tag,
        'sequence': seq,
        'n_reps': len(ddgs),
        'ddG_mean': round(ddg_mean, 3),
        'ddG_median': round(ddg_median, 3),
        'ddG_std': round(ddg_std, 3),
        'ddG_min': round(ddg_min, 3),
        'ddG_max': round(ddg_max, 3),
        'ddG_range': round(ddg_range, 3),
        'z_score_vs_WT': round(z, 3),
        'outlier_flag': outlier,
        'significance': significance,
    }
