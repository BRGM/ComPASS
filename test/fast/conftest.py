import sys
import os


def pytest_configure(config):

    # configure forked option
    try:
        import pytest_forked
    except ImportError:
        raise Exception("missing plugin: pytest-forked")
    config.option.forked = True
