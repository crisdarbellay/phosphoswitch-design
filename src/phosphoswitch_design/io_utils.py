"""
io_utils.py — FASTA, CSV, and PDB file utilities.

All I/O that the pipeline stages share lives here so the logic is not
duplicated across five scripts.
"""

from __future__ import annotations
import csv
import os
import re
from typing import Iterator


# ---------------------------------------------------------------------------
# FASTA
# ---------------------------------------------------------------------------
def parse_mpnn_header(header: str) -> dict[str, str]:
    """Parse a LigandMPNN FASTA header into a metadata dict.

    LigandMPNN writes headers like::

        >PDBNAME, id=42, overall_confidence=0.9347, ligand_confidence=0.8812

    All key=value pairs are extracted.  If the first field has no '=' it
    becomes 'name'.
    """
    meta: dict[str, str] = {}
    for part in header.lstrip(">").split(","):
        part = part.strip()
        if "=" in part:
            k, v = part.split("=", 1)
            meta[k.strip()] = v.strip()
        elif "name" not in meta:
            meta["name"] = part
    return meta


def parse_fasta(path: str, skip_template: bool = True) -> Iterator[tuple[str, str]]:
    """Yield (header_line, sequence) pairs from a FASTA file.

    Parameters
    ----------
    path : str
        Path to the FASTA file.
    skip_template : bool
        If True (default), skip the template entry (id=0 in LigandMPNN output).
        Design sequences have ``id=`` in their header with a non-zero id.
    """
    cur_hdr: str | None = None
    cur_seq = ""

    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if cur_hdr is not None:
                    if not skip_template or "id=" in cur_hdr:
                        yield cur_hdr, cur_seq
                cur_hdr = line
                cur_seq = ""
            elif line:
                cur_seq += line

    if cur_hdr is not None:
        if not skip_template or "id=" in cur_hdr:
            yield cur_hdr, cur_seq


def write_fasta(path: str, records: list[tuple[str, str]]) -> None:
    """Write a list of (header, sequence) pairs to a FASTA file.

    Headers should NOT include the leading '>'; it is added automatically.
    """
    with open(path, "w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n{seq}\n\n")


# ---------------------------------------------------------------------------
# Path metadata extraction (LigandMPNN output tree)
# ---------------------------------------------------------------------------
def get_subspace_from_path(fa_path: str) -> str:
    """Extract subspace name from a path like
    ``.../phase2/<subspace>/T<temp>/<file>.fa``.
    """
    parts = fa_path.replace("\\", "/").split("/")
    try:
        idx = parts.index("phase2")
        return parts[idx + 1]
    except (ValueError, IndexError):
        return "unknown"


def get_state_temp_from_path(fa_path: str) -> tuple[str, str]:
    """Extract MPNN state ('A' or 'B') and temperature string from path."""
    state = "?"
    temp = "?"
    m = re.search(r"lmpnn_out_([AB])", fa_path)
    if m:
        state = m.group(1)
    m = re.search(r"/T([\d.]+)/", fa_path)
    if m:
        temp = m.group(1)
    return state, temp


# ---------------------------------------------------------------------------
# CSV
# ---------------------------------------------------------------------------
def load_csv_indexed(path: str, key: str = "tag") -> dict[str, dict]:
    """Load a CSV file into a dict indexed by the column *key*.

    Rows where the key is empty are silently skipped.
    Missing file returns an empty dict (non-fatal — later stages may be
    partially complete).
    """
    if not os.path.exists(path):
        print(f"  [io_utils] Skipping (not found): {path}")
        return {}
    d: dict[str, dict] = {}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            tag = row.get(key)
            if tag:
                d[tag] = row
    return d


def write_csv(path: str, rows: list[dict], fieldnames: list[str] | None = None) -> None:
    """Write a list of dicts to a CSV file.

    If *fieldnames* is None, keys are taken from the first row (preserving
    insertion order as of Python 3.7+).
    """
    if not rows:
        return
    if fieldnames is None:
        fieldnames = list(rows[0].keys())
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


# ---------------------------------------------------------------------------
# PDB helpers (thin wrappers around sequence_gen functions for convenience)
# ---------------------------------------------------------------------------
def strip_phospho_to_apo(
    phospho_pdb: str,
    output_pdb: str,
    phospho_resid: int = 30,
    phospho_type: str = "TYR",
) -> None:
    """Strip phosphate atoms from a PDB and rename PTR/TYS → TYR.

    Used when an apo backbone must be derived from a phospho PDB, e.g.
    when the only available structure is a phosphorylated form.

    Parameters
    ----------
    phospho_pdb : str
        Input PDB containing the phospho-Tyr (PTR residue or ATOM + HETATM
        phosphate atoms at *phospho_resid*).
    output_pdb : str
        Path to write the apo PDB.
    phospho_resid : int
        Residue number of the phospho-Tyr in the PDB chain (default 30).
    phospho_type : str
        Residue name to use for the de-phosphorylated residue (default TYR).
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
                    continue  # drop phosphate atoms
                new_line = line[:17] + f"{phospho_type:<3}" + line[20:]
                if new_line.startswith("HETATM"):
                    new_line = "ATOM  " + new_line[6:]
                fout.write(new_line)
            else:
                if line.startswith("HETATM"):
                    line = "ATOM  " + line[6:]
                fout.write(line)
