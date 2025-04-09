"""
HyPhy analysis methods package.
"""

from .base import HyPhyMethod
from .busted import BustedMethod
from .relax import RelaxMethod
from .cfel import CfelMethod
from .fel import FelMethod
from .meme import MemeMethod
from .prime import PrimeMethod
from .registry import HyPhyMethodRegistry

__all__ = [
    'HyPhyMethod',
    'BustedMethod',
    'RelaxMethod',
    'CfelMethod',
    'FelMethod',
    'MemeMethod',
    'PrimeMethod',
    'HyPhyMethodRegistry'
]
