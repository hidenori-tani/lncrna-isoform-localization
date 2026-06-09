import json, os, tempfile
from cswitch.keynumbers import set_number, get_number

def test_set_and_get_roundtrip():
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "kn.json")
        set_number("compartment_switcher_fraction", 0.123, path=path)
        assert get_number("compartment_switcher_fraction", path=path) == 0.123
        with open(path) as f:
            assert json.load(f)["compartment_switcher_fraction"] == 0.123

def test_set_preserves_existing_keys():
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "kn.json")
        set_number("a", 1, path=path)
        set_number("b", 2, path=path)
        assert get_number("a", path=path) == 1
        assert get_number("b", path=path) == 2
