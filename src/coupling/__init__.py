"""Step 3 coupling: drive the reference Fortran RT engine with the Step 2 haze.

The Fortran reads its haze from prescribed-haze tables (``haze_data='presc'``);
:mod:`coupling.presc_haze` writes those tables from a haze specification, either
the observational tables themselves (round-trip validation) or the Step 2
microphysics, so the coupled loop needs no Fortran source changes.
"""

from .presc_haze import (
    HazeBand,
    parse_presc,
    write_presc_haze,
    observational_haze,
    microphysics_haze,
)

__all__ = [
    "HazeBand",
    "parse_presc",
    "write_presc_haze",
    "observational_haze",
    "microphysics_haze",
]
