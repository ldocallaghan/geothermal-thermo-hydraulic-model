"""The text report and the JSON export come from the same result."""
import json

import pytest

import site_evaluation as se

pytestmark = pytest.mark.slow


@pytest.mark.parametrize("key", ["soultz", "larderello", "united-downs", "pannonian", "newberry"])
def test_the_report_agrees_with_the_export(evaluated, key, capsys, tmp_path):
    o = evaluated(se.find_site(key))
    se.report(o)
    text = capsys.readouterr().out
    path = tmp_path / "result.json"
    se.export_result(o, path)
    d = json.loads(path.read_text())
    assert d["site"]["name"] in text
    assert f"STATUS: {d['status']['label']}" in text
    for r in d["status"]["reasons"]:
        assert r in text
    p, c = d["programme"], d["cycle"]
    assert f"pipe {p['pipe']}" in text
    assert f"flow {p['flow']:.0f} {p['flow_unit']}" in text
    assert f"Over the cycle: {c['site_verdict']}; deciding case {c['deciding_case']}" in text
    assert f"drilling fluid {c['drill_SG']:.3f} SG (limit {c['drill_SG_hi']:.3f})" in text
    for k in ("drilling_fluid", "trip_fluid"):
        w = d["windows"][k]
        assert (f"{w['fluid'] / 1000:.3f}" if w["open"] else w["note"]) in text
    assert d["tools"]["tool"] in text
    assert d["versions"]["python"]
