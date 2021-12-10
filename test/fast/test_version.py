import ComPASS


def test_version():
    assert ComPASS.__version__ > "4.3"
    assert ComPASS.__version__ >= "4.3"
    assert ComPASS.__version__ > "4.3.4"
    assert ComPASS.__version__ < "99"


if __name__ == "__main__":
    test_version()
