"""
Regression tests using a dict config instead of a config file path.
Mirrors test_regression.py but passes config= rather than config_path=.
"""
import pathlib
import pytest
from scipy.io import readsav
from KN1DPy.kn1d import kn1d
from KN1DPy.kn1d_lite import kn1d_lite
import numpy as np


CONFIG_DICT = {
    "kinetic_h": {
        "mesh_size": 10,
        "ion_rate": "jh",
        "ci_test": True,
        "alpha_cx_test": False,
        "grid_fctr": 0.3,
        "extra_energy_bins_eV": [],
    },
    "kinetic_h2": {
        "mesh_size": 6,
        "grid_fctr": 0.3,
        "extra_energy_bins_eV": [],
        "ci_test": True,
        "alpha_cx_test": False,
    },
    "collisions": {
        "H2_H_EL": True,
        "H2_H2_EL": True,
        "H2_P_EL": True,
        "H2_P_CX": True,
        "H_H_EL": True,
        "H_P_EL": True,
        "H_P_CX": True,
        "SIMPLE_CX": True,
    },
}

INPUT_DIR = pathlib.Path(__file__).parent / "input"

LABELS = [
    "nH", "nH2", "nHP", "TH", "TH2", "THP",
    "GammaxH", "GammaxH2", "Sion", "SideWallH",
    "SH", "SP", "qxH_total", "qxH2_total", "Lyman", "Balmer",
]


@pytest.fixture(scope="module")
def dict_results(run_in_tmp_dir):
    data_file = INPUT_DIR / "cmod_test_in.sav"
    sav = readsav(str(data_file.resolve()))
    Ti = sav["Ti"] * 1e3
    Te = sav["Te"] * 1e3
    n = sav["n"] * 1e20
    return kn1d(
        x=sav["x"],
        xlimiter=sav["xlimiter"],
        xsep=sav["xsep"],
        GaugeH2=sav["GaugeH2"],
        mu=sav["mu"],
        Ti=Ti,
        Te=Te,
        n=n,
        vxi=sav["vxi"],
        LC=sav["LC"],
        PipeDia=sav["PipeDia"],
        max_gen=100,
        compute_errors=1,
        config=CONFIG_DICT,
    )


@pytest.fixture(scope="module")
def file_results(run_in_tmp_dir, config_path):
    data_file = INPUT_DIR / "cmod_test_in.sav"
    sav = readsav(str(data_file.resolve()))
    Ti = sav["Ti"] * 1e3
    Te = sav["Te"] * 1e3
    n = sav["n"] * 1e20
    return kn1d(
        x=sav["x"],
        xlimiter=sav["xlimiter"],
        xsep=sav["xsep"],
        GaugeH2=sav["GaugeH2"],
        mu=sav["mu"],
        Ti=Ti,
        Te=Te,
        n=n,
        vxi=sav["vxi"],
        LC=sav["LC"],
        PipeDia=sav["PipeDia"],
        max_gen=100,
        compute_errors=1,
        config_path=config_path,
    )


@pytest.mark.parametrize("label", LABELS)
def test_kn1d_dict_matches_file(dict_results, file_results, label):
    np.testing.assert_allclose(
        getattr(dict_results, label),
        getattr(file_results, label),
        rtol=1e-10,
        atol=0,
        err_msg=f"Mismatch for {label}",
    )


@pytest.fixture(scope="module")
def lite_dict_result(run_in_tmp_dir):
    data_file = INPUT_DIR / "cmod_test_in.sav"
    sav = readsav(str(data_file.resolve()))
    Ti = sav["Ti"] * 1e3
    Te = sav["Te"] * 1e3
    n = sav["n"] * 1e20
    x = sav["x"]
    xsep = float(sav["xsep"])
    vxi = sav["vxi"]
    mask = x >= xsep
    x_in   = np.concatenate([[0.0],  x[mask] - xsep])
    Ti_in  = np.concatenate([[np.interp(xsep, x, Ti)],  Ti[mask]])
    Te_in  = np.concatenate([[np.interp(xsep, x, Te)],  Te[mask]])
    n_in   = np.concatenate([[np.interp(xsep, x, n)],   n[mask]])
    vxi_in = np.concatenate([[np.interp(xsep, x, vxi)], vxi[mask]])
    return kn1d_lite(
        x=x_in, mu=int(sav["mu"]), Ti=Ti_in, Te=Te_in, n=n_in,
        vxi=vxi_in,
        incident_n0=1e15,
        energies_eV=[3.0],
        config=CONFIG_DICT,
    )


@pytest.fixture(scope="module")
def lite_file_result(run_in_tmp_dir, config_path):
    data_file = INPUT_DIR / "cmod_test_in.sav"
    sav = readsav(str(data_file.resolve()))
    Ti = sav["Ti"] * 1e3
    Te = sav["Te"] * 1e3
    n = sav["n"] * 1e20
    x = sav["x"]
    xsep = float(sav["xsep"])
    vxi = sav["vxi"]
    mask = x >= xsep
    x_in   = np.concatenate([[0.0],  x[mask] - xsep])
    Ti_in  = np.concatenate([[np.interp(xsep, x, Ti)],  Ti[mask]])
    Te_in  = np.concatenate([[np.interp(xsep, x, Te)],  Te[mask]])
    n_in   = np.concatenate([[np.interp(xsep, x, n)],   n[mask]])
    vxi_in = np.concatenate([[np.interp(xsep, x, vxi)], vxi[mask]])
    return kn1d_lite(
        x=x_in, mu=int(sav["mu"]), Ti=Ti_in, Te=Te_in, n=n_in,
        vxi=vxi_in,
        incident_n0=1e15,
        energies_eV=[3.0],
        config_path=config_path,
    )


def test_kn1d_lite_dict_matches_file(lite_dict_result, lite_file_result):
    np.testing.assert_allclose(
        lite_dict_result.nH,
        lite_file_result.nH,
        rtol=1e-10,
        atol=0,
        err_msg="kn1d_lite nH mismatch between dict and file config",
    )
