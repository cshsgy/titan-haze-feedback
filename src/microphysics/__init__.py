"""Titan haze microphysics: fractal-aggregate coagulation + sedimentation
scaling law (Step 2 of the titan-haze-feedback project).

See ``docs/scaling_law.md`` for the derivation.
"""

from .constants import AerosolParams, DEFAULT
from .atmosphere import Atmosphere
from .scaling_law import solve_scaling_law, ScalingResult
from . import transport

__all__ = [
    "AerosolParams",
    "DEFAULT",
    "Atmosphere",
    "solve_scaling_law",
    "ScalingResult",
    "transport",
]
