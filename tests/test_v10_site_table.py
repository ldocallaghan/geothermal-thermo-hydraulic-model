"""The full v1.0 four-site comparative table, frozen (spec test group 1).

Marked slow because it runs Models 1-5 per site and so needs src/water_table.npz
(a one-off ~2 min IAPWS build). With the table cached the four sites take about a
second. Regenerate the golden file with `python tests/make_golden.py`
only when a change to v1.0 behaviour is intended and explained in CHANGELOG.md.

Model 1 quantities (bit temperature, MW, ROP) come out of a BVP solve and move
at the 1e-3 level with numpy/scipy versions, so they are compared loosely.
Verdict strings and the stability flags are compared exactly -- those are what
the README publishes.
"""
import json
import os

import pytest

from conftest import GOLDEN

GOLDEN_FILE = os.path.join(GOLDEN, "v10_site_table.json")

pytestmark = pytest.mark.slow

LOOSE = ("T_bit_C", "MW_prod", "ROP_m_hr", "rop_gain")
EXACT = ("regime", "grc_regime_cold", "breaks", "frac_limited", "survives",
         "stability", "verdict")
TIGHT = ("z_m", "T_rock_C", "anisotropy", "Sv_MPa", "Shmin_MPa", "SHmax_MPa",
         "P_fluid_MPa", "sigma_theta_MPa", "mc_cold_MPa", "P_need_MPa",
         "overbalance_MPa", "m_min_kg_s")


@pytest.fixture(scope="module")
def actual():
    import make_golden
    return {r["name"]: r for r in (make_golden.row(s) for s in make_golden.SITES)}


@pytest.fixture(scope="module")
def expected():
    if not os.path.exists(GOLDEN_FILE):
        pytest.skip("golden file missing; run python tests/make_golden.py")
    with open(GOLDEN_FILE, encoding="utf-8") as f:
        return {r["name"]: r for r in json.load(f)["sites"]}


def test_all_four_sites_present(actual, expected):
    assert set(actual) == set(expected)


def test_v10_table_unchanged(actual, expected):
    problems = []
    for name, exp in expected.items():
        got = actual[name]
        for key in EXACT:
            if got[key] != exp[key]:
                problems.append(f"{name}.{key}: {got[key]!r} != {exp[key]!r}")
        for key in TIGHT:
            if got[key] != pytest.approx(exp[key], rel=1e-9, abs=1e-9):
                problems.append(f"{name}.{key}: {got[key]} != {exp[key]}")
        for key in LOOSE:
            if got[key] != pytest.approx(exp[key], rel=5e-3):
                problems.append(f"{name}.{key}: {got[key]} != {exp[key]}")
    assert not problems, "v1.0 behaviour changed:\n  " + "\n  ".join(problems)


def test_published_verdicts_are_the_readme_verdicts(expected):
    """What the README's comparative table claims, independent of the solver."""
    verdicts = {name.split("(")[0].split("/")[0].strip(): r["verdict"]
                for name, r in expected.items()}
    assert verdicts == {
        "Upper Rhine Graben": "CONDITIONAL",
        "Larderello": "GO",
        "United Downs": "NO-GO (stability)",
        "Pannonian Basin": "GO",
    }
