"""
HyPhy analysis methods package.
"""

from .base import HyPhyMethod
from .busted import BustedMethod
from .cfel import CfelMethod
from .fel import FelMethod
from .meme import MemeMethod
from .prime import PrimeMethod
from .registry import HyPhyMethodRegistry
from .relax import RelaxMethod

__all__ = [
    "HyPhyMethod",
    "BustedMethod",
    "RelaxMethod",
    "CfelMethod",
    "FelMethod",
    "MemeMethod",
    "PrimeMethod",
    "HyPhyMethodRegistry",
]
