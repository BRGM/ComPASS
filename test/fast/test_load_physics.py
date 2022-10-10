import ComPASS


def test_load_diphasic():
    simulation = ComPASS.load_physics("diphasic")
    assert simulation.physics_name() == "diphasic"


def test_load_linear_water():
    simulation = ComPASS.load_physics("linear_water")
    assert simulation.physics_name() == "linear_water"


def test_load_water2ph():
    simulation = ComPASS.load_physics("water2ph")
    assert simulation.physics_name() == "water2ph"
