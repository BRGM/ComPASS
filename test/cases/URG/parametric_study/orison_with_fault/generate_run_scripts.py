import os
import pandas as pd
from os.path import isfile, join


# Reading the parameters file
parameters_file = "parameters.csv"
parameters_array = pd.read_csv(parameters_file, delimiter=";", dtype=str)


# Setting the outputs directory (creating it if not existant)
output_dir = "sensitivity_study_scripts"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)


# Retrieving each parameter value
for i, line in parameters_array.iterrows():
    phi = line["phi"]
    k_socle = line["Ks"]
    k_middle = line["Km"]
    k_top = line["Kt"]
    k_fracture_socle = line["Kfs"]
    k_fracture_meso = line["Kfm"]
    k_fracture_tert_quat = line["Kft"]
    lambda_meso = line["kth_m"]
    lambda_tert_quat = line["kth_t"]
    fracture_thickness = line["ef"]
    # add parameters ... (modify LatinHypercubeSampling.py which provides the parameters.csv)

    print(i)

    # Reading the template file
    templatename = (
        "orison_avec_faille.py"  # Replace the template file with the studied one
    )
    with open(templatename, "r") as file:
        lines = file.readlines()

    # Creating offspring scripts from the template with the different samples of parameters
    output_filename = f"./{output_dir}/sample_{i}.py"
    output_folder_csv = f"/home/haouchine/orison_avec_faille/output-sample_{i}"
    dataframe_cvs = (
        f"/home/haouchine/orison_avec_faille/output-sample_{i}/Data_sample_{i}.csv"
    )
    with open(output_filename, "w") as file1:
        for iline in lines:
            dictionnary = {
                "K_socle": k_socle,
                "K_meso": k_middle,
                "K_tert_quat": k_top,
                "kth_meso": lambda_meso,
                "kth_tert_quat": lambda_tert_quat,
                "K_Fault_socle": k_fracture_socle,
                "K_Fault_meso": k_fracture_meso,
                "K_Fault_tert_quat": k_fracture_tert_quat,
                "f_thickness": fracture_thickness,
                "Num": i,
                "output_folder": output_folder_csv,
                "data_csv": dataframe_cvs,
            }
            mystring_mod = iline
            for keyname, newname in dictionnary.items():
                mystring_mod = mystring_mod.replace(keyname, str(newname))
            file1.write(mystring_mod)


# Running every offspring script
path = r"./sensitivity_study_scripts/"
files = [f for f in os.listdir(path) if isfile(join(path, f))]
print("ok")
# print(files)

num_samples = parameters_array.shape[0]

for i in range(num_samples):
    print(f"Running file sample_{i}.py ...")
    cmd = f"python ./sensitivity_study_scripts/sample_{i}.py"
    os.system(cmd)
