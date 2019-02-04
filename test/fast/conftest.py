def pytest_configure(config):
    try:
        import pytest_forked
    except ImportError:
        raise Exception("missing plugin: pytest-forked")
    config.option.forked = True
