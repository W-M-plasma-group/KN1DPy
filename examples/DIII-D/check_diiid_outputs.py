"""
DIII-D KN1D output checker — Python vs IDL
===========================================
Run from the KN1DPy root directory:

    python examples/DIII-D/check_diiid_outputs.py
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from compare_outputs import run_comparison

run_comparison(
    case_name  = 'DIII-D',
    py_dir     = 'examples/DIII-D/python_output',
    idl_dir    = 'examples/DIII-D/IDL_output',
    idl_prefix = 'diiid_test',
)
