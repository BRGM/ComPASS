

class SimulationContext:
    """
    A structure that holds global parameter.

    It has three attributes:
    - `abort_on_ksp_failure`: will abort the simulation if the iterative linear solver fails, defaults to `False`.
    - `dump_system_on_ksp_failure`: will dump the linear system if the iterative linear solver fails, defaults to `False`.
    - `abort_on_newton_failure`: will abort the simulation if a Newton loop fails (cf. :class:`ComPASS.newton.Newton`), defaults to `False`.
    """

    def __init__(self):
        self.abort_on_ksp_failure = False
        self.dump_system_on_ksp_failure = False
        self.abort_on_newton_failure = False
