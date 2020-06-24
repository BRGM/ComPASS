import sys
import os
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

filename = sys.argv[1]
production_temperatures = np.load(filename)

plt.clf()
times = np.array([data[0] for data in production_temperatures])
temperatures = np.array([data[1] for data in production_temperatures])
plt.plot(times, temperatures - 273.15)  # convert to celsius degrees
plt.xlabel("time (s)")
plt.ylabel("production temperature (Celsius degree)")
plt.savefig(os.path.splitext(filename)[0] + ".png")
