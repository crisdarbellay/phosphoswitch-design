"""
phosphoswitch_design — multi-state protein design of bidirectional phosphoswitches.

Pipeline stages:
    01  LigandMPNN sequence generation across 4 hypothesis tracks
    02  Plausibility filter + geometry-based mechanism scoring
    03  Diversity-capped top-candidate selection with H1/H2 contact filter
    04  Deep PyRosetta 4-state validation (N=20 replicates)
    08  Consensus ranking + wet-lab candidate export

Target: LMNA Y45 phosphorylated by Src kinase; ±29 aa construct (59 aa),
        Y45 = position 30 in the construct.
"""

__version__ = "0.1.0"
__author__ = "Cris Darbellay"

from .mechanism import (
    PHOS_BINDERS,
    PHOS_REPELLERS,
    parse_pdb,
    bidirectional_score,
    mechanism_score,
    classify_h1_h2,
)
from .filtering import (
    passes_plausibility,
    get_mutation_signature,
    passes_h1_h2_filter,
    select_with_diversity,
)
from .io_utils import (
    parse_fasta,
    parse_mpnn_header,
    load_csv_indexed,
)

__all__ = [
    "PHOS_BINDERS",
    "PHOS_REPELLERS",
    "parse_pdb",
    "bidirectional_score",
    "mechanism_score",
    "classify_h1_h2",
    "passes_plausibility",
    "get_mutation_signature",
    "passes_h1_h2_filter",
    "select_with_diversity",
    "parse_fasta",
    "parse_mpnn_header",
    "load_csv_indexed",
]
