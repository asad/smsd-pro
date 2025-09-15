import os, pytest

@pytest.fixture(autouse=True, scope="session")
def _ensure_out_dir():
    OUT = os.path.join(os.path.dirname(os.path.dirname(__file__)), "test_output")
    os.makedirs(OUT, exist_ok=True)
    return OUT
