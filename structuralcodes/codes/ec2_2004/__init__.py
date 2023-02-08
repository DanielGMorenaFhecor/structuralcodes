"""EUROCODE 2 1992-1-1:2004"""
import typing as t

from ._crack_control import (
    crack_min_steel_area,
    crack_min_steel_area_with_prestresed_tendons,
    k_crack_min_steel_area,
    kc_crack_min_steel_area_flanges,
    kc_crack_min_steel_area_pure_tension,
    kc_crack_min_steel_area_rectangular,
    w_max,
)

__all__ = [
    'w_max',
    'crack_min_steel_area',
    'kc_crack_min_steel_area_pure_tension',
    'kc_crack_min_steel_area_rectangular',
    'kc_crack_min_steel_area_flanges',
    'crack_min_steel_area_with_prestresed_tendons',
]

__title__: str = 'EUROCODE 2 1992-1-1'
__year__: str = '2004'
__materials__: t.Tuple[str] = ('concrete',)
