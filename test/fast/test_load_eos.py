import ComPASS


def test_load_diphasic():
    simulation = ComPASS.load_eos("diphasic")
    assert simulation.eos_name() == "diphasic"


def test_load_linear_water():
    simulation = ComPASS.load_eos("linear_water")
    assert simulation.eos_name() == "linear_water"


def test_load_water2ph():
    simulation = ComPASS.load_eos("water2ph")
    assert simulation.eos_name() == "water2ph"
