"""
The worked example for a new site: build the site from the template,
validate it, evaluate it, and write the result as JSON.

    python examples/new_site.py              # validate, and check against Soultz
    python examples/new_site.py --evaluate   # also evaluate (several minutes)
    python examples/new_site.py --evaluate --json soultz.json

The template holds Soultz's inputs, so the site it builds is Soultz's, and
with --evaluate its row of the site table is compared with the pinned one
(tests/golden/v133_site_table.json).
"""
import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "src"))

import comparative_sites as cs        # noqa: E402
import site_evaluation as se          # noqa: E402
from site_template import make_site   # noqa: E402

GOLDEN = os.path.join(HERE, "..", "tests", "golden", "v133_site_table.json")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[1])
    ap.add_argument("--evaluate", action="store_true", help="evaluate the site")
    ap.add_argument("--json", help="write the result to this JSON file")
    a = ap.parse_args()

    site = make_site()
    for w in se.validate_site(site):
        print("warning:", w)
    print(f"{site.name}: inputs valid; {site.tier}")
    print("same inputs as the repository's Soultz:", site == se.SOULTZ)
    if not (a.evaluate or a.json):
        return
    o = se.evaluate(site)
    se.report(o)
    if a.json:
        se.export_result(o, a.json)
        print(f"result written to {a.json}")
    if os.path.exists(GOLDEN):
        pinned = json.load(open(GOLDEN))[se.SOULTZ.name]
        now = cs.golden_rows_v133([o])[site.name]
        same = {k: now[k] == pinned[k] for k in ("status", "verdicts", "site_verdict")}
        print("the Soultz row reproduced:", all(same.values()), same)


if __name__ == "__main__":
    main()
