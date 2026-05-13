from dataclasses import asdict

from KN1DPy.kn1d import kn1d
from scipy.io import readsav
import numpy as np

import argparse
import os
import sys
import time


def run():
    np.set_printoptions(linewidth=225)
    np.set_printoptions(threshold=sys.maxsize)

    standard_out = sys.stdout

    ##Input
    data_file = "./sav_files/kn1d_test_inputs.sav"
    # data_file = './sav_files/1090904018_950to1050.sav'
    # data_file = './sav_files/1090904029_950to1050_towall.sav'
    print("Loading file: " + data_file)
    sav_data = readsav(data_file)

    ##Output
    print("Beginning KN1D")
    start = time.time()
    results = kn1d(
        sav_data["x"],
        sav_data["x_lim"],
        sav_data["x_sep"],
        sav_data["p_wall"],
        sav_data["mu"],
        sav_data["t_i"],
        sav_data["t_e"],
        sav_data["n_e"],
        sav_data["vx"],
        sav_data["lc"],
        sav_data["d_pipe"],
        max_gen=100,
        Hdebug=0,
        H2debug=0,
        debrief=1,
        Hdebrief=1,
        H2debrief=1,
        compute_errors=1,
    )
    end = time.time()
    print("Elapsed Time: ", end - start)
    print()

    output = open("Results/output.txt", "w")
    sys.stdout = output

    for key, value in asdict(results).items():
        print(key)
        print(value)
        print()

    output.close()
    sys.stdout = standard_out


def run_lite():
    """Entry point for running kn1d_lite from an external environment via subprocess."""

    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from KN1DPy.kn1d_lite import kn1d_lite

    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--config", required=True)
    args = parser.parse_args()

    # Load inputs
    d = np.load(args.input)

    x = d["x"]
    mu = float(d["mu"].item())
    Ti = d["Ti"]
    Te = d["Te"]
    n = d["n"]
    vxi = d["vxi"]
    incident_n0 = float(d["incident_n0"].item())
    velocities_ms = d["velocities_ms"]
    fractions = d["fractions"]

    # Run
    res = kn1d_lite(
        x=x,
        mu=mu,
        Ti=Ti,
        Te=Te,
        n=n,
        vxi=vxi,
        incident_n0=incident_n0,
        velocities_ms=velocities_ms,
        fractions=fractions,
        config_path=args.config,
    )

    # Save outputs
    np.savez(
        args.output,
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
