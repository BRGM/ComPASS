import sys
from pathlib import Path
import matplotlib.pyplot as plt
from yaml_utilities import *


def find_divergence_between_logs(dir1, dir2):

    data1 = safe_load(dir1)
    data2 = safe_load(dir2)

    simlog1 = ConvergenceLog(data1)
    simlog2 = ConvergenceLog(data2)
    i = 0
    # simlog.time_step contains a dictionnary
    for (key1, ts1), (key2, ts2) in zip(
        simlog1.time_step.items(), simlog2.time_step.items()
    ):
        breakpoint()
        if ts1 != ts2:
            print(f"Timeline divergence at time step {key1} of the 1st simulation")
            print(f"        and {key2} of the 2nd simulation.\n")
            print(f"time step of 1st simulation:\n{ts1}")
            print(f"time step of 2nd simulation:\n{ts2}")
            exit()
        i += 1


def plot_comparison(directories, x, y):
    """
    Give as argument the pathes to the output directories of the scripts you want to plot.
    The distinct simulations will be plot on the same graph, gathered by convergence
    informations.

    plot_keys(
        'time',
        'newton_iterations_per_attempt',
        'nb_newton_it', 'nb_it_in_success_newton',
        'lsolver_iterations_per_newton',
        'mean_lsolver_iterations_per_newton',
        'success_dt',
        'computation_time')


    """
    plt.clf()
    save = []
    for i, directory in enumerate(directories):
        directory = Path(directory)
        if not directory.is_dir():
            raise IOError(f"Could not find {str(directory)} directory.")

        save.append(directory.stem)
        print(f"Ploting figure number {i} using directory :{directory} ")
        log_dict = analysed_timeloops_log(directory)

        assert x in log_dict, "First argument nor found in dictionary"
        assert y in log_dict, "Second argument nor found in dictionary"

        if x == "newton_iterations_per_attempt" or y == "newton_iterations_per_attempt":
            print("""ERROR: USE OF 'newton_iterations_per_attempt' NOT ALLOWED""")
            return
        plt.plot(log_dict[x], log_dict[y], label="{}".format(directory.stem))
        print()
        print()
        print()
    plt.xlabel(x)
    plt.ylabel(y)
    path = Path().absolute()
    output_dir = path / "comparison"
    output_dir.mkdir(exist_ok=True)
    directories = "-".join(save)
    output_dir = output_dir / f"{directories}"
    output_dir.mkdir(exist_ok=True)

    figfile = output_dir / f"{str(y)}-functions-of-{str(x)}"
    plt.savefig(figfile)
    print(f"Figure {figfile} has been created")


if __name__ == "__main__":
    directories = sys.argv[1:]
    # if(len(directories)==2):
    #     find_divergence_between_logs(sys.argv[1], sys.argv[2])

    plot_comparison(directories, "time", "success_dt")
