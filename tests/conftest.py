"""Shared fixtures. Puts src/ on the path so tests import the models directly."""
import os
import sys

SRC = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

GOLDEN = os.path.join(os.path.dirname(os.path.abspath(__file__)), "golden")
