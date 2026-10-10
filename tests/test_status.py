"""The status of an evaluation: each component's, and the top-level label."""
import copy

import pytest

import site_evaluation as se


def components(**changes):
    c = dict(
        solver=dict(ok=True, converged=True, runs_ok=True, fluid_passes=3),
        tool=dict(ok=True, detail="PASS"),
        coating=dict(status="within rating", pipe="TK-Drakon coated"),
        hole_cleaning=dict(ok=True, v_ann_min=1.6, conflict=None),
        pumps=dict(ok=True, spp=20e6, within_rig=True, gpm=600.0),
        admissibility=dict(inadmissible=[]),
        stability=dict(verdict="GO", states=dict(drilling="GO", connection="GO",
                                                 **{"connection, fresh rock": "GO"}, trip="GO")))
    for k, v in changes.items():
        c[k] = {**c[k], **v}
    return c


def test_all_well_is_go():
    assert se.site_status(components()) == dict(label="GO", reasons=[])


def test_a_solver_failure_is_indeterminate_whatever_else():
    """A run that did not converge leaves every result unknown, even a NO-GO."""
    c = components(solver=dict(ok=False, runs_ok=False),
                   stability=dict(verdict="NO-GO", states=dict(trip="NO-GO")))
    s = se.site_status(c)
    assert s["label"] == "INDETERMINATE"
    assert "did not converge" in s["reasons"][0]
    s = se.site_status(components(solver=dict(ok=False, converged=False)))
    assert s["reasons"] == ["the fluid iteration did not converge"]


def test_a_tool_exceedance_with_no_schedule_is_no_go_with_its_reason():
    c = components(tool=dict(ok=False, detail="FAIL: running in needs continuous circulation"))
    s = se.site_status(c)
    assert s["label"] == "NO-GO"
    assert s["reasons"] == ["tool: running in needs continuous circulation"]


def test_no_go_names_the_failing_states():
    c = components(stability=dict(verdict="NO-GO [CONDITIONAL..NO-GO]",
                                  states=dict(drilling="NO-GO", connection="CONDITIONAL",
                                              trip="NO-GO")))
    s = se.site_status(c)
    assert s["label"] == "NO-GO"
    assert s["reasons"] == ["the wall is not held: drilling, trip"]


def test_a_coating_over_its_rating_is_indeterminate_but_not_over_a_no_go():
    over = dict(coating=dict(status="INDETERMINATE"))
    assert se.site_status(components(**over))["label"] == "INDETERMINATE"
    c = components(**over, tool=dict(ok=False, detail="FAIL: 207 C at the bit while drilling"))
    assert se.site_status(c)["label"] == "NO-GO"


def test_conditional_lists_its_conditions():
    c = components(stability=dict(verdict="CONDITIONAL"), tool=dict(detail="PASS with the schedule found"),
                   pumps=dict(within_rig=False), hole_cleaning=dict(ok=False))
    s = se.site_status(c)
    assert s["label"] == "CONDITIONAL"
    assert s["reasons"] == ["a fluid heavier than water", "the circulation schedule the tool needs",
                            "more pump capacity than the reference rig",
                            "hole cleaning below the FORGE annular velocity"]


@pytest.mark.slow
def test_a_forced_solver_failure_in_an_evaluation_is_indeterminate(evaluated):
    o = copy.copy(evaluated(se.SOULTZ))
    o["m1"] = dict(o["m1"], success=False)
    s = se.site_status(se.component_statuses(o))
    assert s["label"] == "INDETERMINATE"


def test_the_prescription_down_the_hole_is_a_condition():
    """A site GO at the bottom whose common window down the hole needs a fluid
    heavier than water is CONDITIONAL on that fluid; an empty common window
    is a condition too, not a NO-GO."""
    pres = dict(drilling_fluid=dict(open=True, fluid=1036.0, note=None, heavier_than_water=True),
                trip_fluid=dict(open=False, fluid=None, note="none: the scenarios need different fluids",
                                heavier_than_water=False))
    c = components()
    c["prescription"] = pres
    s = se.site_status(c)
    assert s["label"] == "CONDITIONAL"
    assert s["reasons"] == ["a drilling fluid of 1.036 SG at surface",
                            "trip fluid: none: the scenarios need different fluids"]
