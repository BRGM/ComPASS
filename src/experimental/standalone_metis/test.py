import numpy as np
from metis import metis_part_graph as part_graph

connectivity = np.load("connectivity.npz")
neighbors = connectivity["neighbors"]
offsets = connectivity["offsets"]

cells = np.unique(neighbors)
print(neighbors)
print(offsets)
print(cells)
assert np.all(cells == np.arange(cells.shape[0]) + 1)  # Fortran indexing

for k in range(2, min(16, len(cells))):
    print(f"{k} parts - Kway")
    colors = part_graph(neighbors - 1, offsets, k, Kway=True)  # use C-style indexing
    print(np.unique(colors))
    if len(np.unique(colors)) != k:
        print("WARNING", len(np.unique(colors)), "colors vs", k, "parts")

for k in range(2, min(16, len(cells))):
    print(f"{k} parts - recursive")
    colors = part_graph(neighbors - 1, offsets, k)  # use C-style indexing
    print(np.unique(colors))
    assert len(np.unique(colors)) == k
