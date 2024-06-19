import numpy as np
import matplotlib.pyplot as plt
import os
import re
import glob
from scipy.interpolate import RBFInterpolator
import seaborn as sns


# Function to extract number from string
def extract_number(string):
    match = re.search(r"\d+", string)
    return int(match.group()) if match else 0


# Function to perform a Radial Basis Function interpolation
def RBF_interpolation(param1_values, param2_values, misfit):
    # Perform RBF interpolation
    points = np.column_stack((param1_values, param2_values))
    values = np.array(misfit)
    rbf_interpolator = RBFInterpolator(points, values, kernel="cubic", smoothing=0.01)

    # Grid definition for the interpolation
    param1_grid = np.linspace(min(param1_values), max(param1_values), 700)
    param2_grid = np.linspace(min(param2_values), max(param2_values), 700)
    X, Y = np.meshgrid(param1_grid, param2_grid)
    grid_points = np.column_stack((X.ravel(), Y.ravel()))

    # Interpolate misfit values on the grid
    interpolated_misfit = rbf_interpolator(grid_points).reshape(X.shape)

    return X, Y, interpolated_misfit


# Function to ask user which parameter to consider
def ask_input(prompt, valid_range):
    while True:
        try:
            value = int(input(prompt))
            if value in valid_range:
                return value
            else:
                print(f"Please enter a number in the range {valid_range}.")
        except ValueError:
            print("Invalid input. Please enter a number.")


# Load wells data
wells_file_path = "/home/haouchine/orison_avec_faille/Observables_sampling_noise.txt"
Wells_dataframe = np.loadtxt(wells_file_path, skiprows=1, delimiter=" ", dtype=float)

# Extract wells coordinates and temperature
x_coords = Wells_dataframe[:, 0]
y_coords = Wells_dataframe[:, 1]
depth = Wells_dataframe[:, 2]
T_forage = Wells_dataframe[:, 3]

# Initialize useful arrays
mean_T_misfit, Q1_T_misfit, Q3_T_misfit, simulation_numbers = [], [], [], []

# Parent folder
parent_directory = "/home/haouchine/orison_avec_faille/"

# Browse csv files in the parent folder
csv_files = glob.glob(os.path.join(parent_directory, "output-sample_*/*.csv"))
csv_files.sort(key=extract_number)

for csv_file_path in csv_files:

    print("Processing file : ", csv_file_path)

    # Load data frome CSV
    data_table = np.loadtxt(csv_file_path, skiprows=1, delimiter=";", dtype=str)

    # Select the relevant lines
    selected_data = data_table[:15]

    # Compute the misfit
    T_pred = selected_data[:, 7].astype(float)
    T_misfit = np.abs(T_forage - T_pred)

    # Compute useful stats
    moyenne = np.mean(T_misfit)
    Q1 = np.percentile(T_misfit, 25)
    Q3 = np.percentile(T_misfit, 75)

    # Extract the simulation number from the string name
    simulation_number = extract_number(os.path.basename(csv_file_path))

    # Append data
    mean_T_misfit.append(moyenne)
    Q1_T_misfit.append(Q1)
    Q3_T_misfit.append(Q3)
    simulation_numbers.append(simulation_number)


figures_output = "/home/haouchine/orison_avec_faille/figures_output/"
if not os.path.exists(figures_output):
    os.makedirs(figures_output)

# Plotting the samples comparison
plt.figure(figsize=(20, 12))
plt.tick_params(axis="both", which="major", labelsize=6)  # Ajust labels size
plt.scatter(
    simulation_numbers,
    mean_T_misfit,
    color="blue",
    label="Mean Misfit $T_{forage} - T_{simulation}$",
    s=30,
)  # s for the size of the points
plt.scatter(
    simulation_numbers,
    Q1_T_misfit,
    color="green",
    label="Q1 Misfit $T_{forage} - T_{simulation}$",
    s=30,
)
plt.scatter(
    simulation_numbers,
    Q3_T_misfit,
    color="red",
    label="Q3 Misfit $T_{forage} - T_{simulation}$",
    s=30,
)
plt.xlabel("Simulation number")
plt.ylabel("Misfit $T_{forage} - T_{simulation}$")
# plt.xticks(simulation_numbers, rotation='vertical') # Use simulation number as label
# plt.grid(True)
plt.legend(
    framealpha=0.5,
    loc="upper left",
    bbox_to_anchor=(1.05, 1),
    borderaxespad=0,
    fancybox=True,
    shadow=True,
    borderpad=0,
    fontsize=11,
)  # Paramètres de la légende
plt.savefig(figures_output + "Comparaison_series.png", bbox_inches="tight")
plt.close()

# Plot the probability density function of the misfit
plt.figure(figsize=(10, 6))
sns.kdeplot(mean_T_misfit, fill=True)
plt.xlabel("Mean Temperature Misfit")
plt.ylabel("Density")
plt.title("Probability Density Function of Mean Temperature Misfit")
plt.savefig(figures_output + "PDF_Mean_Temperature_Misfit.png", bbox_inches="tight")
plt.close()

# Determine the best (lowest) misfit
min_misfit = np.min(mean_T_misfit)
min_misfit_sample = np.argmin(mean_T_misfit)
print()
print("The min misfit is obtained with the sample n°", min_misfit_sample)
print("The value of the misfit is : ", min_misfit)

# Import the parameters of simulation
parameters_array = np.loadtxt(
    "/home/haouchine/orison_avec_faille/parameters.csv",
    skiprows=1,
    delimiter=";",
    dtype=str,
)
parameter_names = [
    "Phi",
    "Ks",
    "Km",
    "Ktq",
    "Kfs",
    "Kfm",
    "Kftq",
    "kthm",
    "kthtq",
    "ef",
]
Phi = parameters_array[:, 1].astype(float)  # porosity
Ks = parameters_array[:, 2].astype(float)  # Socle permeability
Km = parameters_array[:, 3].astype(float)  # Mesozoic permeability
Ktq = parameters_array[:, 4].astype(float)  # Tertiary-Quaternary permeability
Kfs = parameters_array[:, 5].astype(float)  # Socle fault permeability
Kfm = parameters_array[:, 6].astype(float)  # Mesozoic fault permeability
Kftq = parameters_array[:, 7].astype(float)  # Tertiary-Quaternary fault permeability
kthm = parameters_array[:, 8].astype(float)  # Mesoizoic thermal conductivity
kthtq = parameters_array[:, 9].astype(float)  # Tertiary-Quaternary thermal conductivity
ef = parameters_array[:, 10].astype(float)  # Fault thickness

# Mapping from names to actual values arrays ----------------------------------------------------------------------
parameter_values = {
    "Phi": Phi,
    "Ks": Ks,
    "Km": Km,
    "Ktq": Ktq,
    "Kfs": Kfs,
    "Kfm": Kfm,
    "Kftq": Kftq,
    "kthm": kthm,
    "kthtq": kthtq,
    "ef": ef,
}

# Choose parameters and generate plots
while True:
    print("\nAvailable parameters for study:")
    print()
    print("0. --Exit")

    for i, name in enumerate(parameter_names):
        print(f"{i + 1}. {name}")
    print()

    param1_index = (
        ask_input(
            "Enter the number for the parameter you want to study (or 0 to exit): ",
            range(0, len(parameter_names) + 1),
        )
        - 1
    )

    if param1_index == -1:
        break

    param2_index = (
        ask_input(
            "Enter the number for the second parameter: ",
            range(0, len(parameter_names) + 1),
        )
        - 1
    )

    if param2_index == -1:
        break

    param1 = parameter_names[param1_index]
    param2 = parameter_names[param2_index]

    param1_values = parameter_values[param1]
    param2_values = parameter_values[param2]

    # Perform interpolation
    X_grid, Y_grid, interpolated_misfit = RBF_interpolation(
        param1_values, param2_values, mean_T_misfit
    )

    # Determine if axes should be logarithmic
    log_x = param1 not in ["Phi", "ef", "kthm", "kthtq"]
    log_y = param2 not in ["Phi", "ef", "kthm", "kthtq"]

    # Parameters scatter plot
    plt.figure(figsize=(20, 12))
    scatter = plt.scatter(
        param1_values, param2_values, c=mean_T_misfit, cmap="rainbow", s=30
    )
    plt.colorbar(scatter, label="Mean Temperature Misfit")
    if log_x:
        plt.xscale("log")
    if log_y:
        plt.yscale("log")
    plt.xlabel(param1)
    plt.ylabel(param2)
    plt.title(f"Mean Misfit as a function of {param1} and {param2}")
    plt.savefig(
        f"{figures_output}{param1}_{param2}_comparison.png", bbox_inches="tight"
    )
    plt.close()

    # 3D scatter plot without interpolation
    fig = plt.figure(figsize=(30, 24))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot(
        param1_values,
        param2_values,
        mean_T_misfit,
        linestyle="--",
        color="gray",
        linewidth=1,
    )
    scatter = ax.scatter(
        param1_values,
        param2_values,
        mean_T_misfit,
        c=mean_T_misfit,
        cmap="rainbow",
        s=150,
    )
    ax.set_xlabel(param1)
    ax.set_ylabel(param2)
    ax.set_zlabel("Mean Temperature Misfit")
    ax.set_title("Mean Misfit as a function of " + param1 + " and " + param2)
    fig.colorbar(scatter, label="Mean Temperature Misfit")
    azimut = 20
    élévation = 40
    ax.view_init(azimut, élévation)
    plt.savefig(
        f"{figures_output}{param1}_{param2}_3D_scatter.png", bbox_inches="tight"
    )
    plt.close()

    # 2D Heatmap plot
    plt.figure(figsize=(20, 12))
    contourf = plt.contourf(
        X_grid, Y_grid, interpolated_misfit, cmap="rainbow", levels=100
    )
    contour = plt.contour(
        X_grid, Y_grid, interpolated_misfit, colors="black", levels=10, linewidths=0.5
    )
    plt.colorbar(contourf, label="Mean Temperature Misfit")
    plt.xlabel(param1)
    plt.ylabel(param2)
    if log_x:
        plt.xscale("log")
    if log_y:
        plt.yscale("log")
    plt.title(f"Interpolated Mean Misfit as a function of {param1} and {param2}")
    plt.clabel(contour, inline=True, fontsize=8)
    plt.savefig(
        f"{figures_output}{param1}_{param2}_2D_interpolation_with_contours.png",
        bbox_inches="tight",
    )
    plt.close()

    # 3D plot
    fig = plt.figure(figsize=(20, 12))
    ax = fig.add_subplot(111, projection="3d")
    surf = ax.plot_surface(X_grid, Y_grid, interpolated_misfit, cmap="rainbow")
    ax.set_xlabel(param1)
    ax.set_ylabel(param2)
    ax.set_zlabel("Mean Temperature Misfit")
    ax.set_title(f"Interpolated Mean Misfit as a function of {param1} and {param2}")
    fig.colorbar(surf, label="Mean Temperature Misfit")
    azimut = 20
    élévation = 40
    ax.view_init(azimut, élévation)
    plt.savefig(
        f"{figures_output}{param1}_{param2}_3D_interpolation.png", bbox_inches="tight"
    )
    plt.close()

print("Exiting the program.")
