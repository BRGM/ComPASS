import ComPASS


def test_version():
    assert ComPASS.version > "4.3"
    assert ComPASS.version >= "4.3"
    assert ComPASS.version > "4.3.4"
    assert ComPASS.version < "99"


if __name__ == "__main__":
    test_version()
