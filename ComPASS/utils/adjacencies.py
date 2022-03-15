import numpy as np


def filter_adjacency_table(indices):
    split_points = np.cumsum([len(a) for a in indices])[:-1]
    used, inverse = np.unique(np.hstack(indices), return_inverse=True)
    return used, np.split(inverse, split_points)
