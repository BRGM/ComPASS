import numpy as np
import matplotlib.pyplot as plt
import yaml, sys, glob


"""
This script plots the strong scaling (T(np)/T(1p)) of the computing times spent in
the make_one_timestep function, the Jacobian computation, the matrix filling
and linear solve for different number of processes.
It takes in a directory containing ComPASS output directories named 1p/, 2p/ ... (by the number of procs they were run on)
Each of these directories have to contain a profiling.yaml file with the following structure:

time_step_1:
  total_time_step_time: 201.9858598448336
  jacobian_build_times:
  - 7.5603496711701155
  - 6.587741019204259
  - 6.6055031362921
  matrix_fill_times:
  - 1.641370676457882
  - 0.9255065098404884
  - 0.9325351044535637
  linear_solve_times:
  - 57.366618540138006
  - 51.844311302527785
  - 51.90861861407757
time_step_2:
  total_time_step_time: 109.7446683421731
  jacobian_build_times:
  - 6.599095797166228
  - 6.60384213924408
  - 6.613600347191095
  matrix_fill_times:
  - 0.9186121225357056
  - 0.9223169628530741
  - 0.900626776739955
  linear_solve_times:
  - 2.1487002186477184
  - 51.926561802625656
  - 17.204795315861702
...

Times are computed as the sum of the local time on each proc.
The second argument of the script is the reference (start) case
(Then the scaling plot will be truncated for smaller values and be computed
as T(np)/T(ref))
"""
try:
    basedir = sys.argv[1]
    start_at = float(sys.argv[2])
except IndexError:
    print("Syntax: python3 scaling.py <path/to/output_dirs> <np_start>")
    exit()


def n_procs(path):
    return int(path.split("/")[-1][:-1])


output_dirs = sorted(glob.glob(basedir + "/*"), key=n_procs)
name = " ".join(output_dirs[0].split("/")[-4:-1])
data = []

for path in output_dirs:
    with open(path + "/profiling.yaml") as f:
        if n_procs(path) >= start_at:
            data.append((n_procs(path), yaml.safe_load(f)))


total_times = []
linear_solve_times = []
jacobian_times = []

for n, n_dict in data:
    my_jac_t, my_lin_t, my_tot_t = 0.0, 0.0, 0.0
    for time_dict in n_dict.values():
        # Dividing by nprocs to get the mean time spent by each proc
        my_tot_t += time_dict["total_time_step_time_ext"] / n
        my_jac_t += sum(time_dict["jacobian_build_times"]) / n
        my_lin_t += sum(time_dict["linear_solve_times"]) / n
    total_times.append(my_tot_t)
    jacobian_times.append(my_jac_t)
    linear_solve_times.append(my_lin_t)

tt = total_times[0]
jt = jacobian_times[0]
lst = linear_solve_times[0]
total_speed_up = [tt / time for time in total_times]
jacobian_speed_up = [jt / time for time in jacobian_times]
linear_solve_speed_up = [lst / time for time in linear_solve_times]
xticks = [n_procs(path) for path in output_dirs if n_procs(path) >= start_at]
nc = len(xticks)
theory_reference = [t / xticks[0] for t in xticks]

print("Total times :", total_times)
print("Jacob times :", jacobian_times)
print("Lsolv times :", linear_solve_times)
print(xticks)

fig, (ax, tabax) = plt.subplots(
    nrows=2,
    figsize=(13, 9),
    constrained_layout=True,
    gridspec_kw=dict(height_ratios=[3, 1]),
)
tabax.axis("off")
ax.set_xscale("linear")
ax.set_yscale("linear")
ax.set_xticks(xticks)
# ax.set_xticklabels(xticks)
ax.set_xlabel("Number of procs")
ax.set_ylabel("T(1p)/T(np)")
tab_colors = ["w", "tab:blue", "tab:orange", "tab:olive"]
plot_colors = ["k", "tab:blue", "tab:orange", "tab:olive"]
ax.plot(
    xticks, total_speed_up, label=f"Global speed up", color=plot_colors[0], marker="o"
)
ax.plot(
    xticks,
    jacobian_speed_up,
    label=f"Jacobian speed up",
    color=plot_colors[1],
    marker="o",
)
ax.plot(
    xticks,
    linear_solve_speed_up,
    label=f"Linear solving speed up",
    color=plot_colors[2],
    marker="o",
)
ax.plot(
    xticks,
    theory_reference,
    linestyle="dashed",
    color="grey",
    label=f"Theoretical linear speed up",
)

formatt = lambda times: [f"{time:.1f}" for time in times]
cell_txt = [
    formatt(total_speed_up),
    formatt(jacobian_speed_up),
    formatt(linear_solve_speed_up),
]
the_table = tabax.table(
    cellText=cell_txt,
    cellColours=[[tab_colors[i]] * nc for i in range(3)],
    rowLabels=["Global speed up", "Jacobian speed up", "Linear solving speed up"],
    colLabels=xticks,
    loc="center",
)

ax.set_title(f"Scaling plot : {name}")
the_table.set_fontsize(10)
the_table.scale(1, 1.5)
ax.legend()
plt.savefig("scaling.png")
