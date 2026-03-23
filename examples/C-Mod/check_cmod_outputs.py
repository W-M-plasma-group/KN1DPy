"""
C-Mod KN1D output checker — Python vs IDL
==========================================
Run from the KN1DPy root directory:

    python examples/C-Mod/check_cmod_outputs.py
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from compare_outputs import run_comparison

run_comparison(
    case_name  = 'C-Mod',
    py_dir     = 'examples/C-Mod/python_output',
    idl_dir    = 'examples/C-Mod/IDL_output',
    idl_prefix = 'cmod_test',
)
