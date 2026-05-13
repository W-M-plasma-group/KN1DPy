import pytest


LABELS = [
    "nH",
    "nH2",
    "nHP",
    "TH",
    "TH2",
    "THP",
    "GammaxH",
    "GammaxH2",
    "Sion",
    "SideWallH",
    "SH",
    "SP",
    "qxH_total",
    "qxH2_total",
    "Lyman",
    "Balmer",
]


@pytest.mark.parametrize("results", ["cmod", "diiid", "mastu"], indirect=True)
@pytest.mark.parametrize("label", LABELS)
def test_results(results, label, snaptolshot):
    snaptolshot.assert_allclose(getattr(results, label), rtol=1e-3, atol=1e-3)
