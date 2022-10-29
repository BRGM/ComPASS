import numpy as np
import vtkwriters as vtkw

# open binary file generated in tovtu, contains the complete mesh
mesh = np.load("mesh.npz")
vertices = mesh["vertices"]
hexas = mesh["hexas"]
domain = mesh["domain"]

# function to extract a submesh:
# new numbering of cells and nodes, and save the mapping
def submesh_extraction(cellmask):
    new_vertices_id, new_hexas = np.unique(hexas[cellmask], return_inverse=True)
    new_vertices = vertices[new_vertices_id]
    new_hexas.shape = -1, 8
    return new_vertices, new_hexas, domain[cellmask], new_vertices_id


# identify nodes of a cell domain
def nodes_identification(cellmask):
    selected_nodes_id = np.unique(hexas[cellmask])
    return selected_nodes_id


# identify gallery cells and nodes
# Gallery cells can be removed from the mesh,
# necessity to flag the nodes (before any elimination) to impose CL on them
gallery_domain = [4, 5]  # should be 1, 2, 3, 5, 6, 9 ?
gallery_cells = np.zeros(hexas.shape[0], dtype=bool)
for flag in gallery_domain:
    gallery_cells[domain == flag] = True
gallery_nodes_id = nodes_identification(gallery_cells)
gallery_nodes = np.zeros(len(vertices), dtype=bool)
gallery_nodes[gallery_nodes_id] = True

# compute cell centers
def compute_cell_centers(vertices, hexas):
    centers = np.copy(vertices[hexas[:, 0]])
    for j in range(1, 8):
        centers += vertices[hexas[:, j]]
    centers *= 1 / 8
    return centers


centers = compute_cell_centers(vertices, hexas)
x = centers[:, 0]
y = centers[:, 1]
# new selection depending on cell center coordinates
# to extract a small mesh
selection = (x > 665) & (x < 678) & (y > 250) & (y < 312)
print(f"{selection.sum()} cells have been selected")

# this will overwrite global variables
vertices, hexas, domain, new_vertices_id = submesh_extraction(selection)
centers = centers[selection]
gallery_nodes = gallery_nodes[new_vertices_id]

# compute shape ratio
n = hexas.shape[0]
# shape ratio = ||center - vertex|| / min(half edge length)
# cf. https://clas.ucdenver.edu/math-clinic/sites/default/files/attached-files/tech-x_final_report.pdf p.7
half_edge_length = np.empty((n, 3), dtype="d")
shape_ratio = np.empty(n, dtype="d")

for k, hexa in enumerate(hexas):
    v = vertices[hexa]
    half_edge_length[k, 0] = 0.5 * abs(v[1, 0] - v[0, 0])
    half_edge_length[k, 1] = 0.5 * abs(v[2, 1] - v[0, 1])
    half_edge_length[k, 2] = 0.5 * abs(v[4, 2] - v[0, 2])
    shape_ratio[k] = np.linalg.norm(v[0] - centers[0]) / min(half_edge_length[k])


def write_vtu(vertices, hexas, domain, shape_ratio, filename):
    vtkw.write_vtu(
        vtkw.vtu_doc(
            vertices,
            hexas,
            celldata={"domain": domain, "shape ratio": shape_ratio},
        ),
        filename,
    )


# write_vtu(vertices, hexas, domain, shape_ratio, filename = "original_small_mesh.vtu")

remove_shrink_cells = True
if remove_shrink_cells:
    # remove cells which have half_edge_length < threshold
    threshold = 0.1
    print(f"Processing {n} cells with edge threshold {threshold}")
    print("Minimum edges:", np.min(half_edge_length, axis=0))
    print("Maximum shape ratio:", np.max(shape_ratio))

    keep = np.ones(n, dtype=bool)
    # new indices
    ni = np.arange(vertices.shape[0], dtype="i")

    def shrink_along(axis, edges):
        dropped = half_edge_length[:, axis] < threshold
        keep[dropped] = False
        for hexa in hexas[dropped]:
            for i, j in edges:
                vertices[ni[hexa[i]]] += vertices[ni[hexa[j]]]
                vertices[ni[hexa[i]]] *= 0.5
                ni[hexa[j]] = ni[hexa[i]]
        axis_name = {0: "Ox", 1: "Oy", 2: "Oz"}[axis]
        print("Removing", np.sum(dropped), f"cells along {axis_name}")

    edges = [(0, 1), (2, 3), (4, 5), (6, 7)]
    shrink_along(0, edges)
    edges = [(0, 3), (1, 2), (4, 7), (5, 6)]
    shrink_along(1, edges)
    edges = [(0, 4), (1, 5), (2, 6), (3, 7)]
    shrink_along(2, edges)

    kept_vertices, inverse = np.unique(ni, return_inverse=True)
    vertices = vertices[kept_vertices]
    hexas = inverse[ni[hexas[keep]]]
    domain = domain[keep]
    shape_ratio = shape_ratio[keep]
    gallery_nodes = gallery_nodes[kept_vertices]

    print(f"{hexas.shape[0]} cells kept after removing shrink cells")

gallery_cells = np.zeros(hexas.shape[0], dtype=bool)
for flag in gallery_domain:
    gallery_cells[domain == flag] = True
# new selection to remove the gallery cells
selection = np.invert(gallery_cells)
if gallery_cells.any():
    with_gal_n_nodes = vertices.shape[0]
    nogal_vertices, nogal_hexas, nogal_domain, new_vertices_id = submesh_extraction(
        selection
    )
    nogal_shape_ratio = shape_ratio[selection]
    gallery_nodes = gallery_nodes[new_vertices_id]

gallery_nodes_id = np.flatnonzero(gallery_nodes)
cell_mapping = np.zeros(selection.size, dtype=int)
cell_mapping[selection] = range(nogal_hexas.shape[0])
node_mapping = np.zeros(with_gal_n_nodes, dtype=int)
node_mapping[new_vertices_id] = range(nogal_vertices.shape[0])

np.savez(
    "small_mesh_with_gal.npz",
    vertices=vertices,
    hexas=hexas,
    domain=domain,
    cell_mapping=cell_mapping,
    node_mapping=node_mapping,
)
write_vtu(
    vertices,
    hexas,
    domain,
    shape_ratio,
    filename="small_mesh_with_gal.vtu",
)

print(f"{nogal_hexas.shape[0]} cells kept after removing gallery cells")
np.savez(
    "small_mesh_no_gal.npz",
    vertices=nogal_vertices,
    hexas=nogal_hexas,
    domain=nogal_domain,
    gallery_nodes_id=gallery_nodes_id,
)
write_vtu(
    nogal_vertices,
    nogal_hexas,
    nogal_domain,
    nogal_shape_ratio,
    filename="small_mesh_no_gal.vtu",
)
