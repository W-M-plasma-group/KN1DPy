"""
Entry point for running kn1d_lite from an external environment via subprocess.

Usage:
    /path/to/KN1DPy/.pixi/envs/default/bin/python run_kn1d_lite.py \
        --input /path/to/input.npz --output /path/to/output.npz --config /path/to/config.toml

See docs/kn1d_lite.md for full usage instructions.
"""

import argparse
import json
import os
import sys
import tempfile
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from KN1DPy.kn1d_lite import kn1d_lite

parser = argparse.ArgumentParser()
parser.add_argument('--input',  required=True)
parser.add_argument('--output', required=True)
parser.add_argument('--config', required=True)
args = parser.parse_args()

# Load inputs
d = np.load(args.input)

x           = d['x']
mu          = float(d['mu'].item())
Ti          = d['Ti']
Te          = d['Te']
n           = d['n']
vxi         = d['vxi']
incident_n0 = float(d['incident_n0'].item())
velocities_ms = d['velocities_ms']
fractions     = d['fractions']

# Run
res = kn1d_lite(
    x=x, mu=mu, Ti=Ti, Te=Te, n=n, vxi=vxi,
    incident_n0=incident_n0,
    velocities_ms=velocities_ms,
    fractions=fractions,
    config_path=args.config,
)

# Save outputs
np.savez(args.output,
         xH=res.xH,
         fH=res.fH,
         nH=res.nH,
         GammaxH=res.GammaxH,
         VxH=res.VxH,
         TH=res.TH,
         qxH_total=res.qxH_total,
         Sion=res.Sion,
         fHBC=res.fHBC,
         GammaxHBC=np.array([res.GammaxHBC]),
         vr=res.vr,
         vx=res.vx,
         Tnorm=np.array([res.Tnorm]),
)
