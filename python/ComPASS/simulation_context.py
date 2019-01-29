import ComPASS

class SimulationContext:

    def __init__(self):
        self.abort_on_ksp_failure = False
        self.dump_system_on_ksp_failure = False
        self.abort_on_newton_failure = False


