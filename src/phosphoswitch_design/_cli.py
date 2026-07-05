"""
_cli.py — entry points for console_scripts.

Each function simply delegates to the corresponding script's main().
This allows the pipeline stages to be invoked either as:
    python scripts/01_generate.py ...
or, after pip install, as:
    psw-generate ...
"""

import sys
from pathlib import Path

# Ensure the scripts directory is on the path when called via entry point
_SCRIPTS_DIR = Path(__file__).resolve().parents[3] / "scripts"


def _run_script(script_name: str) -> None:
    """Import and run a pipeline script's main() function."""
    import importlib.util
    script_path = _SCRIPTS_DIR / script_name
    if not script_path.exists():
        sys.exit(f"Script not found: {script_path}")
    spec = importlib.util.spec_from_file_location("_script", script_path)
    mod = importlib.util.module_from_spec(spec)  # type: ignore[arg-type]
    spec.loader.exec_module(mod)  # type: ignore[union-attr]
    mod.main()


def generate() -> None:
    """psw-generate — stage 01: LigandMPNN sequence generation."""
    _run_script("01_generate.py")


def filter_seqs() -> None:
    """psw-filter — stage 02: plausibility filter + mechanism scoring."""
    _run_script("02_filter.py")


def select() -> None:
    """psw-select — stage 03: top diverse candidate selection."""
    _run_script("03_select.py")


def rosetta() -> None:
    """psw-rosetta — stage 04: deep 4-state Rosetta validation."""
    _run_script("04_rosetta.py")


def final() -> None:
    """psw-final — stage 08: consensus ranking + wet-lab selection."""
    _run_script("08_select_final.py")
