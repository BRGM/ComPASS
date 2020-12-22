from ComPASS.timestep_management import FixedTimeStep
from ComPASS.timestep_management import TimeStepManager


def test_fixed_timestep():
    print(list(FixedTimeStep(1).steps()))


def test_classical_timestep():
    tsm = TimeStepManager(1, decrease_factor=0.1)
    print(list(tsm.steps()))
    tsm.current_step = 2e-6
    print(list(tsm.steps()))


def test_karma_timestep():
    def generate_timesteps(tsm, kth):
        dt = []

        def success():
            for t in tsm.steps():
                dt.append(t)
                break  # timestep *success*

        def failure():
            i = 0
            for t in tsm.steps():
                i += 1
                if i > 2:  # timestep *failure*
                    dt.append(t)
                    break

        success()  # first timestep
        for _ in range(kth + 1):
            success()
        failure()
        for _ in range(kth):
            success()
        return dt

    kth = 3
    tsm = TimeStepManager(
        1, decrease_factor=0.5, increase_factor=2, karma_threshold=kth
    )
    dt = generate_timesteps(tsm, kth)
    print(dt)
    result = [1, 2, 2, 2, 4, 1, 1, 1, 2]
    assert all([int(dtk) == i for dtk, i in zip(dt, result)])
    tsm = TimeStepManager(
        1,
        decrease_factor=0.5,
        increase_factor=2,
        karma_threshold=kth,
        karma_reset=False,
    )
    dt = generate_timesteps(tsm, kth)
    print(dt)
    result = [1, 2, 4, 8, 16, 8, 8, 8, 16]
    assert all([int(dtk) == i for dtk, i in zip(dt, result)])
