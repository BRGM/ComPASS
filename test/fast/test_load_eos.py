import ComPASS


def test_load_diphasic():
    ComPASS.load_eos("diphasic")


def test_load_linear_water():
    ComPASS.load_eos("linear_water")


def test_load_water2ph():
    ComPASS.load_eos("water2ph")
