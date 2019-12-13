import sys
import os


def pytest_configure(config):
    # add package directory to sys.path
    conftest_dir = os.path.abspath(os.path.dirname(__file__))
    pkg_dir = os.path.join(conftest_dir, '../..')
    sys.path.insert(0, pkg_dir)

    # configure forked option
    try:
        import pytest_forked
    except ImportError:
        raise Exception("missing plugin: pytest-forked")
    config.option.forked = True
