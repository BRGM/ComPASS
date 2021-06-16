import subprocess as sp
import itertools as it


def test_combinations(script, options, repeat=2):
    for comb in it.combinations(options, repeat):
        args = ["python3", script]
        for opt in comb:
            args += [opt, options[opt]]
        print("\n")
        print("#" * 70)
        print(f"Running ComPASS with command : \n{args}")
        print("#" * 70)
        call = sp.run(args, check=True)
        print(call)


def test_product(script, options1, options2):
    for comb in it.product(options1, options2, repeat=1):
        args = ["python3", script]
        args += [comb[0], options1[comb[0]]]
        args += [comb[1], options2[comb[1]]]
        print("\n")
        print("#" * 70)
        print(f"Running ComPASS with command : \n{args}")
        print("#" * 70)
        call = sp.run(args, check=True)
        print(call)


if __name__ == "__main__":
    test_product(
        "vertical_fracture_from_options.py",
        {
            "--lsolver.new.iterative.pc.bjacobi": "True",
            "--lsolver.legacy.iterative.activate_cpramg": "False",
        },
        {
            "--callbacks.dump_system_on_linear_failure": "True",
            "--callbacks.abort_on_linear_failure": "True",
            "--callbacks.abort_on_newton_failure": "True",
            "--callbacks.abort_on_newton_failure": "True",
            "--lsolver.view": "True",
        },
    )
    test_product(
        "doublet_from_options.py",
        {"--lsolver.new": "True", "--lsolver.legacy": "True"},
        {
            "--callbacks.linear_system_dump": "0.0",
            "--callbacks.linear_system_binary_dump": "0.0",
            "--callbacks.abort": "0.0",
            "--callbacks.newton_log": "newton_log.txt",
        },
    )
    test_combinations(
        "doublet_from_options.py",
        {
            "--lsolver.new.iterative.pc.bjacobi": "True",
            "--lsolver.new.iterative.bcgs": "True",
            "--lsolver.new.iterative.tolerance": "1.e-8",
            "--lsolver.new.iterative.gmres.restart": "50",
            "--lsolver.new.iterative.maxit": "200",
        },
        repeat=4,
    )
    test_combinations(
        "doublet_from_options.py",
        {
            "--lsolver.legacy.iterative.activate_cpramg": "False",
            "--lsolver.legacy.iterative.tolerance": "1.e-8",
            "--lsolver.legacy.iterative.gmres.restart": "50",
            "--lsolver.legacy.iterative.maxit": "200",
        },
        repeat=3,
    )
