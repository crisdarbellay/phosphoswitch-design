"""
sequence_gen.py — LigandMPNN wrapper and PDB pre-processing.

LigandMPNN (v_32_010_25) is called as a subprocess because it does not
expose a stable Python API.  This module handles:

    1. PDB pre-processing
       - split_phospho_to_ligand: extracts phosphate atoms into a P4X HETATM
         record so LigandMPNN treats them as a ligand (required for the
         ligand_mpnn model type)
       - strip_phospho_to_apo: removes phosphate atoms and renames PTR→TYR

    2. Checkpoint discovery
       find_protein_mpnn_checkpoint: locates proteinmpnn_v_48_*.pt for
       no-ligand (apo) runs; falls back to ligand_mpnn checkpoint

    3. Core runner
       run_ligandmpnn: constructs and executes the CLI command, returns
       (success, message)

    4. Resume detection
       output_already_done: checks whether a FASTA with design sequences
       already exists in an output directory

Track definitions (4 hypothesis tracks)
-----------------------------------------
Each track specifies:
    id         short label ('1A', '1B', '2A', '2B')
    label      descriptive name used in output directory names
    hypothesis 'H1' or 'H2' — which backbone phospho stabilises
    phos_pdb   backbone PDB for the PHOSPHO state (used as MPNN target
               for state A runs, with P4X ligand inserted)
    apo_pdb    backbone PDB for the APO state (used for state B runs)
    apo_needs_strip  True if apo_pdb must be derived by stripping a
               phospho PDB (track 1A only)

Hypothesis explanation:
    H1: phospho stabilises STRAIGHT helix
        phos target = straight+phos (stateA_phospho.pdb)
        apo  target = hairpin       (stateA_phospho_pulled.pdb, stripped)
    H2: phospho stabilises HAIRPIN
        phos target = pulled hairpin+phos (stateA_phospho_pulled.pdb)
        apo  target = straight            (stateA_aln.pdb)

Tracks 2A/2B add the B-form natural hairpin as the apo backbone,
testing whether designs that evolved on the authentic B-form are
better-behaved than those designed on the artificially pulled A-form.
"""

from __future__ import annotations
import os
import shutil
import subprocess
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# LigandMPNN parameters — set via config or overridden by CLI
# ---------------------------------------------------------------------------
CHAIN = "A"
PHOSPHO_RESID = 30
PHOSPHO_TYPE = "TYR"
LIGAND_CUTOFF = 14.0   # Å radius for ligand context
BATCH_SIZE = 10

# 4 hypothesis tracks
TRACKS = [
    {
        "id": "1A",
        "label": "StraightPhos_PulledApo",
        "hypothesis": "H1",
        "phos_pdb": "output/phase1/stateA_phospho.pdb",
        "apo_pdb":  "output/phase1/stateA_phospho_pulled.pdb",
        "apo_needs_strip": True,
    },
    {
        "id": "1B",
        "label": "PulledPhos_StraightApo",
        "hypothesis": "H2",
        "phos_pdb": "output/phase1/stateA_phospho_pulled.pdb",
        "apo_pdb":  "output/phase1/stateA_aln.pdb",
        "apo_needs_strip": False,
    },
    {
        "id": "2A",
        "label": "StraightPhos_BHairpinApo",
        "hypothesis": "H1",
        "phos_pdb": "output/phase1/stateA_phospho.pdb",
        "apo_pdb":  "output/phase1/stateB_aln.pdb",
        "apo_needs_strip": False,
    },
    {
        "id": "2B",
        "label": "BHairpinPhos_StraightApo",
        "hypothesis": "H2",
        "phos_pdb": "output/phase1/stateB_phospho.pdb",
        "apo_pdb":  "output/phase1/stateA_aln.pdb",
        "apo_needs_strip": False,
    },
]


# ---------------------------------------------------------------------------
# PDB pre-processing
# ---------------------------------------------------------------------------
def split_phospho_to_ligand(
    phospho_pdb: str,
    output_pdb: str,
    phospho_resid: int = PHOSPHO_RESID,
    phospho_type: str = PHOSPHO_TYPE,
) -> None:
    """Extract phosphate atoms from a phospho-Tyr into a P4X ligand HETATM.

    LigandMPNN requires the phosphate to be a separate ligand record (HETATM
    with residue name P4X in chain X) rather than embedded in the protein
    ATOM records.  Protein atoms are written in residue-number order.

    Atom name mapping:
        P  → P,  O1P/OP1 → O1,  O2P/OP2 → O2,  O3P/OP3 → O3

    Parameters
    ----------
    phospho_pdb : str
        Input PDB with phospho-Tyr at *phospho_resid*.
    output_pdb : str
        Path to write the processed PDB.
    phospho_resid : int
        Residue number of the phospho-Tyr in chain A (default 30).
    phospho_type : str
        Residue name for the de-phosphorylated Tyr ATOM lines (default TYR).
    """
    protein_lines_by_res: dict[int, list[str]] = {}
    phospho_atom_lines: list[str] = []

    with open(phospho_pdb) as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            try:
                resnum = int(line[22:26].strip())
            except ValueError:
                continue
            atom_name = line[12:16].strip()

            if resnum == phospho_resid:
                if atom_name in ('P', 'O1P', 'O2P', 'O3P', 'OP1', 'OP2', 'OP3'):
                    phospho_atom_lines.append(line)
                else:
                    new_line = line[:17] + f"{phospho_type:<3}" + line[20:]
                    if new_line.startswith("HETATM"):
                        new_line = "ATOM  " + new_line[6:]
                    protein_lines_by_res.setdefault(resnum, []).append(new_line)
            else:
                if line.startswith("HETATM"):
                    line = "ATOM  " + line[6:]
                protein_lines_by_res.setdefault(resnum, []).append(line)

    # Build P4X ligand records
    _atom_map = {
        'P': 'P', 'O1P': 'O1', 'O2P': 'O2', 'O3P': 'O3',
        'OP1': 'O1', 'OP2': 'O2', 'OP3': 'O3',
    }
    hetatm_lines: list[str] = []
    for line in phospho_atom_lines:
        atom_name = line[12:16].strip()
        canonical = _atom_map.get(atom_name, atom_name)
        new_name = f" {canonical:<3}"
        new_line = (
            "HETATM" + line[6:12] + new_name + line[16:17]
            + "P4X" + " " + "X" + "   1" + line[26:]
        )
        hetatm_lines.append(new_line)

    with open(output_pdb, "w") as fh:
        for resnum in sorted(protein_lines_by_res.keys()):
            for line in protein_lines_by_res[resnum]:
                fh.write(line)
        fh.write("TER\n")
        for line in hetatm_lines:
            fh.write(line)
        fh.write("END\n")


def strip_phospho_to_apo(
    phospho_pdb: str,
    output_pdb: str,
    phospho_resid: int = PHOSPHO_RESID,
    phospho_type: str = PHOSPHO_TYPE,
) -> None:
    """Strip phosphate atoms and rename PTR/TYS → TYR.

    Used when the apo backbone must be derived from a phospho PDB
    (Track 1A: the pulled hairpin exists only in phospho form).

    Parameters
    ----------
    phospho_pdb : str
        Input PDB containing the phospho-Tyr.
    output_pdb : str
        Path to write the apo PDB.
    phospho_resid : int
        Residue number of the phospho-Tyr (default 30).
    phospho_type : str
        Residue name for the de-phosphorylated residue (default TYR).
    """
    with open(phospho_pdb) as fin, open(output_pdb, "w") as fout:
        for line in fin:
            if not line.startswith(("ATOM", "HETATM")):
                fout.write(line)
                continue
            try:
                resnum = int(line[22:26].strip())
            except ValueError:
                fout.write(line)
                continue
            atom_name = line[12:16].strip()

            if resnum == phospho_resid:
                if atom_name in ('P', 'O1P', 'O2P', 'O3P', 'OP1', 'OP2', 'OP3'):
                    continue
                new_line = line[:17] + f"{phospho_type:<3}" + line[20:]
                if new_line.startswith("HETATM"):
                    new_line = "ATOM  " + new_line[6:]
                fout.write(new_line)
            else:
                if line.startswith("HETATM"):
                    line = "ATOM  " + line[6:]
                fout.write(line)


# ---------------------------------------------------------------------------
# Checkpoint discovery
# ---------------------------------------------------------------------------
def find_protein_mpnn_checkpoint(ligandmpnn_dir: str) -> str | None:
    """Return the path to a proteinmpnn_v_48_*.pt checkpoint, or None.

    Used for apo (no-ligand) design runs.  If no proteinmpnn checkpoint is
    found, the caller falls back to the ligand_mpnn checkpoint.

    Parameters
    ----------
    ligandmpnn_dir : str
        Root directory of the LigandMPNN installation.
    """
    params_dir = Path(ligandmpnn_dir) / "model_params"
    candidates = list(params_dir.glob("proteinmpnn_v_48_*.pt"))
    if candidates:
        return str(candidates[0])
    return None


# ---------------------------------------------------------------------------
# Resume detection
# ---------------------------------------------------------------------------
def output_already_done(out_dir: str) -> bool:
    """Return True if LigandMPNN has already produced FASTAs in *out_dir*.

    Checks for at least one FASTA in ``<out_dir>/seqs/`` that contains a
    design sequence (header with ``id=``).

    Parameters
    ----------
    out_dir : str
        The temperature-level output directory (e.g. ``.../T0.1/``).
    """
    seqs_dir = Path(out_dir) / "seqs"
    if not seqs_dir.is_dir():
        return False
    fastas = list(seqs_dir.glob("*.fa"))
    if not fastas:
        return False
    for fa in fastas:
        with open(fa) as fh:
            n_designs = sum(
                1 for line in fh
                if line.startswith(">") and "id=" in line
            )
        if n_designs > 0:
            return True
    return False


# ---------------------------------------------------------------------------
# Core LigandMPNN runner
# ---------------------------------------------------------------------------
def run_ligandmpnn(
    pdb_path: str,
    out_dir: str,
    has_ligand: bool,
    design_positions: list[int],
    fixed_residues: list[int],
    num_seqs: int,
    temperature: float,
    ligandmpnn_dir: str,
    ligandmpnn_ckpt: str,
    chain: str = CHAIN,
    ligand_cutoff: float = LIGAND_CUTOFF,
    batch_size: int = BATCH_SIZE,
    seed: int = 42,
    omit_aa: list[str] | None = None,
) -> tuple[bool, str]:
    """Run LigandMPNN's run.py on a single (PDB, temperature) combination.

    Parameters
    ----------
    pdb_path : str
        Input PDB.  Must contain P4X HETATM if has_ligand=True.
    out_dir : str
        Output directory.  FASTAs will appear in ``<out_dir>/seqs/``.
    has_ligand : bool
        True → use ligand_mpnn model with P4X ligand context.
        False → use protein_mpnn model (apo backbone runs).
    design_positions : list[int]
        1-indexed residue positions that are designable.
    fixed_residues : list[int]
        1-indexed residue positions that are held fixed (complement of design_positions).
    num_seqs : int
        Total sequences to generate; split into ceil(num_seqs/batch_size) batches.
    temperature : float
        Sampling temperature.  Higher → more diverse / higher mutation rate.
    ligandmpnn_dir : str
        Path to LigandMPNN installation (contains run.py).
    ligandmpnn_ckpt : str
        Path to ligandmpnn_v_32_010_25.pt checkpoint.
    seed : int
        Random seed passed to LigandMPNN (reproducibility).
    omit_aa : list[str] or None
        Amino acid letters to exclude from design (e.g. ['P'] to prevent
        proline from kinking the helix).

    Returns
    -------
    (success, log_message) : (bool, str)
    """
    pdb_abs = os.path.abspath(pdb_path)
    out_abs = os.path.abspath(out_dir)
    os.makedirs(out_abs, exist_ok=True)

    if has_ligand:
        model_type = "ligand_mpnn"
        ckpt_flag = "--checkpoint_ligand_mpnn"
        ckpt = ligandmpnn_ckpt
    else:
        pmpnn_ckpt = find_protein_mpnn_checkpoint(ligandmpnn_dir)
        if pmpnn_ckpt:
            model_type = "protein_mpnn"
            ckpt_flag = "--checkpoint_protein_mpnn"
            ckpt = pmpnn_ckpt
        else:
            model_type = "ligand_mpnn"
            ckpt_flag = "--checkpoint_ligand_mpnn"
            ckpt = ligandmpnn_ckpt

    fixed_str = " ".join(f"{chain}{r}" for r in fixed_residues)
    redesigned_str = " ".join(f"{chain}{r}" for r in design_positions)
    num_batches = max(1, num_seqs // batch_size)

    cmd = [
        sys.executable, os.path.join(ligandmpnn_dir, "run.py"),
        "--model_type", model_type,
        ckpt_flag, ckpt,
        "--pdb_path", pdb_abs,
        "--out_folder", out_abs,
        "--redesigned_residues", redesigned_str,
        "--fixed_residues", fixed_str,
        "--number_of_batches", str(num_batches),
        "--batch_size", str(batch_size),
        "--temperature", str(temperature),
        "--ligand_mpnn_cutoff_for_score", str(ligand_cutoff),
        "--seed", str(seed),
    ]
    if omit_aa:
        cmd += ["--omit_AA", "".join(omit_aa)]

    try:
        result = subprocess.run(
            cmd,
            cwd=ligandmpnn_dir,
            capture_output=True,
            text=True,
            timeout=3600,
        )
        if result.returncode != 0:
            err_tail = "\n".join((result.stderr or "").splitlines()[-15:])
            return False, f"rc={result.returncode}\n{err_tail}"
        return True, "ok"
    except subprocess.TimeoutExpired:
        return False, "timeout (>1h)"
    except Exception as exc:
        return False, f"{type(exc).__name__}: {exc}"


# ---------------------------------------------------------------------------
# Track preparation
# ---------------------------------------------------------------------------
def prepare_track_pdbs(
    track: dict,
    phase2_root: Path,
    phospho_resid: int = PHOSPHO_RESID,
    phospho_type: str = PHOSPHO_TYPE,
) -> tuple[Path, Path]:
    """Pre-process input PDBs for a single hypothesis track.

    Returns
    -------
    (phos_with_ligand, apo_clean) : (Path, Path)
        phos_with_ligand — straight/hairpin+phos PDB with P4X HETATM ligand
        apo_clean        — apo backbone PDB (copy or stripped)
    """
    phos_with_ligand = phase2_root / "phospho_with_ligand.pdb"
    apo_clean = phase2_root / "apo_clean.pdb"

    split_phospho_to_ligand(
        track["phos_pdb"], str(phos_with_ligand), phospho_resid, phospho_type
    )

    if track.get("apo_needs_strip"):
        strip_phospho_to_apo(
            track["apo_pdb"], str(apo_clean), phospho_resid, phospho_type
        )
    else:
        shutil.copy(track["apo_pdb"], apo_clean)

    return phos_with_ligand, apo_clean
