class SimulationInfo:
    def __init__(self):
        self.activate_cpramg = True
        self.activate_direct_solver = False
        self.system = None
        self.ghosts_synchronizer = None
        self.petsc = None


initialized = False
mesh_is_local = False
info = SimulationInfo()
