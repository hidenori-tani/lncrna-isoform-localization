import json, os

DEFAULT_PATH = os.path.join(os.path.dirname(os.path.dirname(__file__)), "key_numbers.json")

def _load(path):
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    return {}

def set_number(key, value, path=DEFAULT_PATH):
    data = _load(path)
    data[key] = value
    with open(path, "w") as f:
        json.dump(data, f, indent=2, sort_keys=True, ensure_ascii=False)

def get_number(key, path=DEFAULT_PATH):
    return _load(path)[key]
