"""The NWG 55-29 cross-check: the stress model rebuilt from Davatzes and
Hickman (2011) and the breakout pattern predicted by rock type, pinned."""
import numpy as np
import pytest

import crosscheck_newberry as nb

PSI = nb.PSI


def test_stresses_at_8420_ft_match_the_paper():
    """Sv by construction; pore pressure from the static survey within 20 psi;
    Shmin at mu 0.70 within 20 psi of the 4,675 psi the paper quotes (its
    formula gives 5,118 at mu 0.55, the value the text calls the best estimate)."""
    z = float(nb.tvd_ft(nb.MD_8420)) * nb.FT
    assert nb.profile(0.55).Sv(z) / PSI == pytest.approx(8852.9, abs=1)
    assert float(nb.profile(0.55).Pp(z)) / PSI == pytest.approx(3121.6, abs=20)
    assert nb.profile(0.70).Shmin(z) / PSI == pytest.approx(4675.3, abs=20)


def test_ucs_porosity_relation():
    assert nb.ucs_from_porosity(0.0) / PSI == pytest.approx(13800.0)
    assert nb.ucs_from_porosity(10.0) / PSI == pytest.approx(13800.0 * np.exp(-0.4744))


@pytest.mark.parametrize("mu", [0.55, 0.70])
def test_pattern_by_rock_type(mu):
    """The volcanic width (modal 36 deg) is matched at 25-35 K of wall
    cooling, but the logged granodiorite breaks out too, where none were
    seen: the pattern is not reproduced."""
    r = nb.crosscheck(mu)
    assert r["consistent_cooling"] == []
    assert r["volcanics_cooling"] == [25.0, 30.0, 35.0]
    for row in r["rows"]:
        if row["cooling"] in r["volcanics_cooling"]:
            assert row["s"]["granodiorite (logged)"]["frac"] > 0.5
