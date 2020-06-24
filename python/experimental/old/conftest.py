import sys
import os

PKG_RELATIVE_PATH = "."


def add_package_to_path():
    conftest_path = os.path.dirname(__file__)
    pkg_path = os.path.abspath(os.path.join(conftest_path, PKG_RELATIVE_PATH))
    sys.path.insert(0, pkg_path)


add_package_to_path()
