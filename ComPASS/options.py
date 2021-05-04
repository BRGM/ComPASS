import sys


def get(name, default=None):
    _, *args = sys.argv
    if name not in args:
        return default
    idx = args.index(name) + 1
    assert idx < len(args), "no value provided after option"
    return args[idx]


def get_bool(name, default=False):
    _, *args = sys.argv
    if name not in args:
        return default
    else:
        return True


if __name__ == "__main__":
    print("toto=", get("--toto"))
    print("tutu=", get("--tutu", "pas de tutu"))
    print("titi=", get("--titi", "pas de titi"))
