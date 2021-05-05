import subprocess as sp

executable = "python3"
script = "solving_from_options.py"
opt_dict = {
    "--dump_ls": ("0.0",),
    "--dump_ls_binary": ("0.0",),
    "--linear_solver_view": ("",),
    "--cpr_amg_type": ("gamg",),
    "--disable_cpramg": ("",),
    "--direct_linear_solver": ("",),
    "--newton_log": ("newton_log.txt",),
}
n = 0
for version in ("legacy", "new"):
    base_args = [executable, script, "--linear_solver_version", version]
    for opt in opt_dict:
        args = [opt] + list(opt_dict[opt])
        print("\n")
        print("#" * 70)
        print(f"Running ComPASS with command : \n{base_args + args}")
        print("#" * 70)
        call = sp.run(base_args + args, check=True)
        print(call)

script = "ksp_failure_vertical_fracture.py"
failure_opts = [
    "--abort_on_linear_failure",
    "--abort_on_newton_failure",
    "--dump_system_on_linear_failure",
]
base_args = [executable, script]
for opt in failure_opts:
    args = base_args + [opt]
    print("\n")
    print("#" * 70)
    print(f"Running ComPASS with command : \n{base_args + args}")
    print("#" * 70)
    call = sp.run(base_args + args, check=True)
    print(call)
