import numpy as np
import matplotlib.pyplot as plt
import yaml, sys, glob

""" This script plots a bar diagram of the relative computing times spent in
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
"""

try:
    basedir = sys.argv[1]
except IndexError:
    print(
        "Please input a directory name containing outputs for different numbers of procs"
    )
    exit()


def n_procs(path):
    return int(path.split("/")[-1][:-1])


output_dirs = sorted(glob.glob(basedir + "/*"), key=n_procs)
name = " ".join(output_dirs[0].split("/")[-4:-1])
data = []

for path in output_dirs:
    with open(path + "/profiling.yaml") as f:
        data.append((n_procs(path), yaml.safe_load(f)))


nc = len(output_dirs)
total_times = []
linear_solve_times = []
jacobian_times = []
matrix_fill_times = []
remainder = []

for n, n_dict in data:
    my_jac_t, my_lin_t, my_tot_t, my_fill_t, my_rem_time = 0.0, 0.0, 0.0, 0.0, 0.0
    for time_dict in n_dict.values():
        # Dividing by nprocs to get the mean time spent by each proc
        my_tot_t += time_dict["total_time_step_time"] / n
        my_jac_t += sum(time_dict["jacobian_build_times"]) / n
        my_lin_t += sum(time_dict["linear_solve_times"]) / n
        my_fill_t += sum(time_dict["matrix_fill_times"]) / n
    my_rem_time += my_tot_t - my_jac_t - my_lin_t - my_fill_t
    total_times.append(my_tot_t)
    jacobian_times.append(my_jac_t)
    linear_solve_times.append(my_lin_t)
    matrix_fill_times.append(my_fill_t)
    remainder.append(my_rem_time)

fig = plt.figure(figsize=(13, 6))
ax = fig.add_subplot(1, 1, 1)
np_ticks = np.arange(len(output_dirs))
np_labels = [f"{n_procs(path)} p" for path in output_dirs]
cell_colors = [
    ["w"] * nc,
    ["tab:blue"] * nc,
    ["tab:red"] * nc,
    ["tab:orange"] * nc,
    ["tab:olive"] * nc,
]

time_array = [
    total_times,
    jacobian_times,
    matrix_fill_times,
    linear_solve_times,
    remainder,
]
proportion = [
    [(jt / tt) * 100 for jt, tt in zip(this, total_times)] for this in time_array[1:]
]
y_offset = np.zeros(nc)
for i, p in enumerate(proportion):
    ax.bar(np_ticks, p, color=cell_colors[i + 1], bottom=y_offset, tick_label="")
    y_offset += p
    for i, v in enumerate(p):
        ax.text(
            i,
            y_offset[i] - 0.5 * v,
            f"{v:.2f}%",
            ha="center",
            va="center_baseline",
            color="k",
        )

formatt = lambda times: [f"{time:.1f}" for time in times]
cell_txt = [
    formatt(total_times),
    formatt(jacobian_times),
    formatt(matrix_fill_times),
    formatt(linear_solve_times),
    formatt(remainder),
]
the_table = ax.table(
    cellText=cell_txt,
    cellColours=cell_colors,
    rowLabels=[
        "Total time (s)",
        "Jacobian computing time (s)",
        "Matrix fill time (s)",
        "Linear solve time (s)",
        "Remaining time (s)",
    ],
    bbox=[0.04, -0.45, 0.95, 0.4],
    colLabels=np_labels,
    loc="bottom",
)

the_table.set_fontsize(10)
the_table.scale(1, 2)
ax.set_ylabel("CPU time repartition (%)")
plt.subplots_adjust(left=0.17, bottom=0.3, right=0.99)
ax.set_title(name)
plt.savefig("bar_cpu_times.png")
