"""Shared fixtures. Puts src/ on the path so tests import the models directly."""
import os
import sys

import pytest

SRC = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

GOLDEN = os.path.join(os.path.dirname(os.path.abspath(__file__)), "golden")


@pytest.fixture(scope="session")
def evaluated():
    """site_evaluation.evaluate, run once per site and set of arguments for the
    whole session. The results are only read by the tests."""
    import site_evaluation as se
    cache = {}

    def get(site, **kw):
        key = (site.name, tuple(sorted(kw.items())))
        if key not in cache:
            cache[key] = se.evaluate(site, **kw)
        return cache[key]

    return get
