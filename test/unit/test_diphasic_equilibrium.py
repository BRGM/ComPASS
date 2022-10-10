import numpy as np

import ComPASS
from ComPASS.utils.units import *
from ComPASS.petrophysics.models.Beaude2018 import Pc


simulation = ComPASS.load_physics("diphasic")
pg = 1 * bar
T = degC2K(20)
Sg = np.linspace(0, 1)

Ca = np.vstack([simulation.diphasic_equilibrium((pg, pg - Pc(Sgk)), T) for Sgk in Sg])


try:
    import matplotlib
except ImportError:
    print("WARNING - matplotlib was not found - no graphics will be generated")
else:
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.clf()
    fig, axes = plt.subplots(1, 3, constrained_layout=True)

    fig.suptitle(f"equilibrium at $p^g$={pg/bar} bar, $T$={K2degC(T)} °C")

    axes[0].plot(Sg, (pg - Pc(Sg)) / bar)
    axes[0].set_xlabel("Sg")
    axes[0].set_ylabel("pl (bar)")

    axes[1].plot(Sg, Ca[:, 0])
    axes[1].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    axes[1].set_xlabel("Sg")
    axes[1].set_ylabel("$C^g_a$")

    axes[2].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    axes[2].plot(Sg, Ca[:, 1])
    axes[2].set_xlabel("Sg")
    axes[2].set_ylabel("$C^l_a$")

    plt.savefig("C_at_equilibrium.png")
