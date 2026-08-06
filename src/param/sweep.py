# -*- coding: utf-8 -*-
"""
Created on Thu Jan 29 17:13:36 2026

@author: k2473520
"""

# v8_split_8/param/sweep.py
from __future__ import annotations
from itertools import product
from typing import Dict, Iterable, List, Mapping


def grid_sweep(base: Mapping[str, float], grid: Mapping[str, Iterable[float]]) -> List[Dict[str, float]]:
    keys = list(grid.keys())
    values = [list(grid[k]) for k in keys]
    out: List[Dict[str, float]] = []
    for combo in product(*values):
        d = dict(base)
        for k, v in zip(keys, combo):
            d[k] = float(v)
        out.append(d)
    return out
