import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.stats import qmc
import pandas as pd


""" Function to generate Quasi Monte Carlo Latin Hypercubes samples for given parameters """


def LH_sampling(n_samples, n_parameters, ranges):

    lhs = qmc.LatinHypercube(
        d=n_parameters, scramble=False
    )  # Setting the Latin Hypercube Sampler (d = dimension)
    samples = lhs.random(n_samples)  # Creating n samples

    for i in range(n_parameters):

        min, max = ranges[i]
        threashold = 1e-10

        # verify the amplitude of the scale
        if max - min < threashold:
            samples[:, i] = 10 ** (
                np.log10(min) + samples[:, i] * (np.log10(max) - np.log10(min))
            )  # logarithmic scale if the amplitude is small
        else:
            samples[:, i] = (
                samples[:, i] * (max - min) + min
            )  # Linear scale if the amplitude is regular

    return samples


""" Ranges of tested parameters """
phi_max = 1
phi_min = 0
Ks_max = 1e-13
Ks_min = 1e-14
Km_max = 1e-12
Km_min = 1e-13
Kt_max = 1e-17
Kt_min = 1e-18
Kfs_max = 1e-12
Kfs_min = 1e-13
Kfm_max = 1e-13
Kfm_min = 1e-14
Kft_max = 1e-15
Kft_min = 1e-16
kth_m_max = 3
kth_m_min = 2
kth_t_max = 2
kth_t_min = 1
f_max = 80
f_min = 70

""" Sampling """
n_samples, n_parameters = 700, 10
ranges = [
    (phi_min, phi_max),
    (Ks_min, Ks_max),
    (Km_min, Km_max),
    (Kt_min, Kt_max),
    (Kfs_min, Kfs_max),
    (Kfm_min, Kfm_max),
    (Kft_min, Kft_max),
    (kth_m_min, kth_m_max),
    (kth_t_min, kth_t_max),
    (f_min, f_max),
]


samples = []
LHS = LH_sampling(n_samples, n_parameters, ranges)

# print("Latin Hypercube Samples:")
for i, sample in enumerate(LHS):
    samples.append(sample)
    # print(f"Sample {i+1}: {sample}")

""" create a dataframe with the sampled parameters values """
parameters = pd.DataFrame(
    samples,
    columns=["phi", "Ks", "Km", "Kt", "Kfs", "Kfm", "Kft", "kth_m", "kth_t", "ef"],
)
print(parameters)

parameters.to_csv("parameters.csv", sep=";")
