\
import os, pytest

def pytest_configure(config):
    out = os.path.join(os.path.dirname(os.path.dirname(__file__)), "test_output")
    os.makedirs(out, exist_ok=True)
