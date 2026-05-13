import os
import pathlib

import pytest

from scipy.io import readsav

from KN1DPy.kn1d import kn1d


@pytest.fixture(scope="session")
def input_directory():
    return pathlib.Path(__file__).parent / "input"


@pytest.fixture(scope="session")
def config_path(input_directory: pathlib.Path):
    return input_directory / "config.toml"


@pytest.fixture(scope="session")
def run_in_tmp_dir(tmp_path_factory):
    """
    Temporarily change to the working directory.
    """
    old_cwd = pathlib.Path.cwd()

    tmp_path = tmp_path_factory.mktemp("kn1d")

    os.chdir(tmp_path)

    yield tmp_path

    os.chdir(old_cwd)


@pytest.fixture(scope="session")
def results(request, run_in_tmp_dir, config_path, input_directory):
    data_file = input_directory / f"{request.param}_test_in.sav"

    sav = readsav(str(data_file.resolve()))

    # Unit conversions (sav file stores keV and 1e20 m^-3)
    Ti = sav["Ti"] * 1e3  # keV  -> eV
    Te = sav["Te"] * 1e3  # keV  -> eV
    n = sav["n"] * 1e20  # 1e20 m^-3 -> m^-3

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
