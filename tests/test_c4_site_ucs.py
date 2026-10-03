"""C4: per-site rock strength (spec v1.1).

sigma_cm takes UCS as an argument, defaulting to the global C.UCS; every site
run passes its own SiteProfile.UCS. The temperature derating is unchanged.
Done when a Cornwall run uses 180 MPa (or the C8-calibrated value, once it
exists). The v1.0 adapter (evaluate(..., v10=True)) keeps the global 200 MPa.
"""
import pytest

import geo_constants as C
import model5_convergence_confinement as m5
import site_evaluation as se
from comparative_sites import CORNWALL, SITES

MPa = 1.0e6


# ------------------------------------------------------------------ sigma_cm
def test_default_is_the_global_ucs():
    for T in (25.0, 180.0, 200.0, 400.0):
        assert m5.sigma_cm(T) == m5.sigma_cm(T, C.UCS)


def test_default_is_read_at_call_time(monkeypatch):
    """None, not a bound default, so C.UCS stays the one source of the global."""
    monkeypatch.setattr(C, "UCS", 150.0 * MPa)
    assert m5.sigma_cm(25.0) == 150.0 * MPa


@pytest.mark.parametrize("T_C, factor", [
    (25.0, 1.0),
    (200.0, 1 - 0.0009 * 175),
    (400.0, 1 - 0.0009 * 375),
    (1000.0, 0.4),              # the floor
])
def test_derating_is_unchanged_and_scales_with_ucs(T_C, factor):
    assert m5.sigma_cm(T_C, 180.0 * MPa) == pytest.approx(180.0 * MPa * factor)
    assert m5.sigma_cm(T_C, 180.0 * MPa) / m5.sigma_cm(T_C) == pytest.approx(0.9)


def test_grc_uses_the_ucs_it_is_given():
    p_i, p0 = 120.0 * MPa, 170.0 * MPa
    assert m5.grc(p_i, p0, 200.0) == m5.grc(p_i, p0, 200.0, UCS=C.UCS)
    weak = m5.grc(p_i, p0, 200.0, UCS=100.0 * MPa)
    strong = m5.grc(p_i, p0, 200.0, UCS=300.0 * MPa)
    assert weak[2] > strong[2]          # weaker rock needs more support (p_cr)


def test_breakout_uses_the_ucs_it_is_given():
    z = CORNWALL.target_depth
    _, Shmin, SHmax, P_fluid = se.stresses_v10(CORNWALL, z)
    T_rock = float(CORNWALL.geotherm(z))
    b = se.breakout_v10(Shmin, SHmax, P_fluid, T_rock, UCS=CORNWALL.UCS)
    assert b["mc_cold"] == pytest.approx(
        m5.KMC * P_fluid + m5.sigma_cm(200.0, 180.0 * MPa))
    # 20 MPa less strength at 25 C, derated to 200 C, over (1 + KMC)
    b10 = se.breakout_v10(Shmin, SHmax, P_fluid, T_rock)
    assert (b["P_need"] - b10["P_need"]) / MPa == pytest.approx(
        20.0 * (1 - 0.0009 * 175) / (1 + m5.KMC), abs=1e-6)


# ------------------------------------------------- full site runs (slow)
@pytest.mark.slow
@pytest.mark.parametrize("site", SITES, ids=lambda s: s.name.split()[0])
def test_every_site_run_passes_its_own_ucs(site):
    """Each site's own strength cases (C8 gives United Downs four; the rest
    have one, their UCS); the row's case carries its UCS through."""
    o = se.evaluate(site)
    assert o["UCS"] in [u for _, u in site.strength_cases]
    assert o["UCS"] == o["stress"]["UCS"]
    assert {c["UCS"] for c in o["stress_cases"]} == {u for _, u in site.strength_cases}
    b = o["breakout"]
    # effective-stress limit since C3: Pp + sigma_cm + KMC (Pw - Pp), with the
    # wall at Model 1's circulating bottom-hole temperature since C7
    Pp, T_wall = o["Pp"], o["T_wall"]
    assert T_wall == pytest.approx(o["m1"]["T_bottom_delivered"])
    assert b["mc_cold"] == pytest.approx(
        Pp + m5.sigma_cm(T_wall, o["UCS"]) + m5.KMC * (o["P_fluid"] - Pp))
    assert o["grc"]["reg_c"] == m5.grc(o["P_fluid"], o["Shmin"], T_wall,
                                       UCS=o["UCS"])[3]


@pytest.mark.slow
def test_cornwall_run_uses_the_c8_calibration_and_the_v10_adapter_does_not():
    """C4's done-when: 180 MPa, 'or the calibrated value from C8 once that
    exists'. It does; the v1.0 adapter keeps the global 200 MPa."""
    o = se.evaluate(CORNWALL)
    assert {c["UCS"] for c in o["stress_cases"]} == {203e6, 172e6, 147e6, 118e6}
    assert se.evaluate(CORNWALL, v10=True)["UCS"] == C.UCS == 200.0 * MPa
