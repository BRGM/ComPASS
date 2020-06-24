import sys
import yaml
import pickle

import numpy as np

from matplotlib import pyplot as plt

with open("analysis.yaml") as f:
    simulations = yaml.load(f, Loader=yaml.FullLoader)
print(f"Loaded {len(simulations)} simulations info.")

for simulation in simulations:
    hash_code = simulation["hash"]
    with open(f"results-{hash_code}", "rb") as f:
        simulation["results"] = pickle.load(f)

markers = ["kx", "r+", "bx", "gx"]
line_markers = ["r-", "g:", "b-", "k:"]
marker = lambda k: markers[k % len(markers)]
line_marker = lambda k: line_markers[k % len(line_markers)]

plt.clf()
for k, simulation in enumerate(simulations):
    params = simulation["params"]
    results = simulation["results"]
    nl = params["nb_layers"]
    plt.plot(results["timesteps"], marker(k), label=f"{nl} layers")
    plt.yscale("log")
plt.legend()
plt.savefig("timesteps.pdf")

plt.clf()
for k, simulation in enumerate(simulations):
    params = simulation["params"]
    results = simulation["results"]
    nl = params["nb_layers"]
    plt.plot(results["times"], results["pwh"], line_marker(k), label=f"{nl} layers")
plt.xscale("log")
plt.legend()
plt.savefig("well head pressure.pdf")

plt.clf()
for k, simulation in enumerate(simulations):
    params = simulation["params"]
    results = simulation["results"]
    nl = params["nb_layers"]
    plt.plot(results["times"], results["qwh"], marker(k), label=f"{nl} layers")
plt.xscale("log")
plt.legend()
plt.savefig("well head flow rate.pdf")

1 / 0

plt.clf()
plt.grid(True)
for k, simulation in enumerate(simulations):
    params = simulation["params"]
    results = simulation["results"]
    x = results["axis_vertices"][:, 0]
    pend = results["axis_pressure"][-1]
    theta = params["theta"]
    plt.plot(x, pend, markers[k], label=f"theta=pi/{int(np.pi/theta)}")
plt.legend()
plt.savefig("final_axis_pressures.pdf")

plt.clf()
plt.grid(True)
markers = ["kx", "r+", "bx", "gx"]
for k, simulation in enumerate(simulations):
    params = simulation["params"]
    results = simulation["results"]
    p = np.array(results["axis_pressure"])
    t = np.array(results["times"])
    vertices = results["axis_vertices"]
    iO = np.argmin(vertices[:, 0])
    pO = p[:, iO]
    theta = params["theta"]
    plt.plot(
        t[t < 10000], pO[t < 10000], markers[k], label=f"theta=pi/{int(np.pi/theta)}"
    )
plt.legend()
plt.savefig("pO.pdf")
