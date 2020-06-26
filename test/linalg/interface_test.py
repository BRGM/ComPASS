from doublet_on_cartesian_grid_setup import *

"""
This script checks important functions of the LinearSolver interface
"""
dt = 86400
kernel = get_kernel()


def test_interface(lsolver):
    newton = Newton(simulation, 1e-5, 8, lsolver)
    newton.init_iteration()
    kernel.Residu_reset_history()
    kernel.Residu_compute(dt)
    newton.convergence_scheme.reset_references(dt)
    kernel.Jacobian_ComputeJacSm(dt)

    print(" * Fill matrix and RHS :")
    lsolver.linear_system.set_from_jacobian()
    print(" * Solve :")
    x, nit, status = lsolver.solve()
    print(f"Nombre d'itérations : {nit}, status : {status}")
    print(" * Residual norm :")
    lsolver.linear_system.check_residual_norm()
    print(" * Ascii dump :")
    lsolver.linear_system.dump_ascii("dump_test/")
    print(" * Binary dump : ")
    lsolver.linear_system.dump_binary("dump_test/")


all_solvers = {
    "Legacy Direct": {"legacy": True, "direct": True},
    "Modern Iterative": {"legacy": False, "direct": False},
    "Modern Direct": {"legacy": False, "direct": True},
    "Error check": {"legacy": False, "direct": True, "activate_cpramg": True},
}

for name, kwargs in all_solvers.items():
    print("\n****************************")
    print("Testing %s Linear Solver\n" % name)
    test_interface(linear_solver(simulation, **kwargs))
