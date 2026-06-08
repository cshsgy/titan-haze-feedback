"""Step 1: 1-D DISORT radiative energy balance for Titan.

Two-band (shortwave solar + longwave thermal) plane-parallel radiative transfer
via pydisort, with haze optical properties supplied by the Step 2 microphysics
and a gray-gas placeholder for the CH4 + CIA opacity (to be replaced with
correlated-k).  Solves for the radiative-equilibrium temperature profile.

Requires the pydisort runtime (see .rtenv); import lazily so the microphysics
package stays importable without torch/pydisort installed.
"""

from .column import Column
from .optics import OpticsParams

__all__ = ["Column", "OpticsParams"]
