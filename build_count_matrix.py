#!/usr/bin/env python3
"""Compatibility wrapper for tools/build_count_matrix.py."""

import runpy
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "tools" / "build_count_matrix.py"
runpy.run_path(str(SCRIPT), run_name="__main__")
