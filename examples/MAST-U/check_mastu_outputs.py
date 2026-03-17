"""
MAST-U KN1D output checker — Python vs IDL
===========================================
Run from the KN1DPy root directory:

    python examples/MAST-U/check_mastu_outputs.py
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from compare_outputs import run_comparison

run_comparison(
    case_name  = 'MAST-U',
    py_dir     = 'examples/MAST-U/python_output',
    idl_dir    = 'examples/MAST-U/IDL_output',
    idl_prefix = 'mastu_test',
)
