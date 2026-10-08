"""Every site input has a stated data basis and source; sites without
measured or calibrated stress and strength, checked against a well, are
classed speculative."""
import pytest

import site_evaluation as se
from comparative_sites import CORNWALL, LARDERELLO, PANNONIAN, SITES

INPUTS = ("stress", "strength", "pore pressure", "temperature", "well check")


@pytest.mark.parametrize("site", SITES, ids=lambda s: s.name.split()[0])
def test_every_input_has_a_basis_and_a_source(site):
    for k in INPUTS:
        basis, source = site.data_basis[k]
        assert basis in se.BASES and source


def test_tiers():
    assert {s.name.split()[0]: s.tier for s in SITES} == {
        "Upper": "evidence-based", "United": "evidence-based",
        "Larderello": "speculative", "Pannonian": "speculative", "Newberry": "speculative"}


def test_newberry_is_held_back_only_by_its_well_check():
    """Its log-derived strength counts; the cross-check not reproducing the
    breakout pattern is what keeps it speculative."""
    from comparative_sites import NEWBERRY
    assert dict(NEWBERRY.missing()) == {"well check": "not reproduced"}
    assert NEWBERRY.mu == 0.55
    assert {p.mu for p in NEWBERRY.stress_cases} == {0.55, 0.70}


@pytest.mark.parametrize("site", [LARDERELLO, PANNONIAN], ids=["Larderello", "Pannonian"])
def test_speculative_sites_name_what_is_missing(site):
    assert dict(site.missing()) == {"stress": "regime bounds", "strength": "unsourced",
                                    "well check": "none"}


def test_soultz_strength_is_the_published_lab_range():
    assert dict(se.SOULTZ.strength_cases) == {"lab UCS 100": 100e6, "lab UCS 115": 115e6,
                                              "lab UCS 130": 130e6}


def test_larderello_follows_venelle_2():
    T = LARDERELLO.temperature()
    assert float(T(2200.0)) == pytest.approx(350.0, abs=1.0)
    assert float(T(2815.0)) == pytest.approx(504.0, abs=1.0)
    assert 507.0 <= float(T(2900.0)) <= 517.0
    assert float(T(LARDERELLO.target())) == pytest.approx(400.0, abs=0.5)


def test_united_downs_is_180_c_at_5_km():
    assert float(CORNWALL.temperature()(5000.0)) == pytest.approx(180.0)
    assert float(CORNWALL.temperature()(CORNWALL.target())) == pytest.approx(400.0)


@pytest.mark.parametrize("site", SITES, ids=lambda s: s.name.split()[0])
def test_v10_adapter_keeps_the_v10_temperatures(site):
    assert site.temperature(v10=True) is site.geotherm
    assert site.target(v10=True) == site.target_depth
