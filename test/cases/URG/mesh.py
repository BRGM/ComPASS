# from collections import namedtuple
import numpy as np
import MeshTools as MT
import vtkwriters as vtkw
import MeshTools.URG as urg

Point = urg.Point
Vector = urg.Vector
Plane = urg.Plane
Hull = urg.Hull


def extract_data(filename):
    data = np.loadtxt(filename, skiprows=1, dtype=np.double)
    index = np.asarray(data[:, 0], dtype=np.intc)
    return data[:, 1:], index


origins, oi = extract_data("URG_origins.out")
normals, ni = extract_data("URG_nvecs.out")

assert np.all(ni == oi)

xmin, ymin, zmin = origins.min(axis=0)
xmax, ymax, zmax = origins.max(axis=0)

Lh, Lv = 1e4, 1e3

# the hull is defined by the intersection of the negative side of an arbitrary
# set of planes
hull = Hull(
    [
        Plane(Point(xmin - Lh, ymin, zmin), Vector(-1, 0, 0)),
        Plane(Point(xmax + Lh, ymin, zmin), Vector(1, 0, 0)),
        Plane(Point(xmin, ymin - Lh, zmin), Vector(0, -1, 0)),
        Plane(Point(xmin, ymax + Lh, zmin), Vector(0, 1, 0)),
        Plane(Point(xmin, ymin, zmin - Lv), Vector(0, 0, -1)),
        Plane(Point(xmin, ymin, zmax + Lv), Vector(0, 0, 1)),
    ]
)

# n = 5
# origins = origins[:n]
# normals = normals[:n]

# output normals as vtu
vtkw.write_vtu(
    vtkw.points_as_vtu_doc(
        origins,
        pointdata={"normals": normals},
    ),
    "normals.vtu",
)

origins = [Point(*a) for a in origins]
normals = [Vector(*a) for a in normals]

mesh = urg.Mesh(hull, list(zip(origins, normals)))

print("built", len(list(mesh.surfaces())), "surface meshes")


def dump_model(model, filename):
    mesh_arrays = [S.as_arrays() for S in model.surfaces()]
    nodes_offset = np.cumsum(
        [
            0,
        ]
        + [a[0].shape[0] for a in mesh_arrays]
    )
    faces_offset = np.cumsum(
        [
            0,
        ]
        + [a[1].shape[0] for a in mesh_arrays]
    )
    all_vertices = np.vstack([a[0] for a in mesh_arrays])
    all_faces = np.vstack(
        [a[1] + offset for a, offset in zip(mesh_arrays, nodes_offset[:-1])]
    )
    fault_id = np.zeros(all_faces.shape[0])
    for fi, t in enumerate(zip(faces_offset[:-1], faces_offset[1:])):
        fault_id[t[0] : t[1]] = fi
    vtkw.write_vtu(
        vtkw.vtu_doc(all_vertices, all_faces, celldata={"id": fault_id}), filename
    )


dump_model(mesh, "all_faults.vtu")

# for k, S in enumerate(mesh.surfaces()):
#    mesh = MT.TSurf.make(*S.as_arrays())


# for k, S in enumerate(mesh.surfaces()):
#    mesh = MT.TSurf.make(*S.as_arrays())
#    MT.to_vtu(mesh, "surface%03d.vtu" % k)

print("adding dtm")

filename = "raster2750.txt"

import Raster

info = Raster.RasterInfo(filename)
raw_data = np.loadtxt(filename, skiprows=6)
data = np.ma.array(raw_data, mask=raw_data == float(info.nodata))
centers = info.centers()

xyz = np.transpose(
    np.vstack(
        [a.ravel() for a in centers]
        + [
            raw_data.ravel(),
        ]
    )
)

xyz = xyz[raw_data.ravel() != float(info.nodata)]

# center dtm
xyz[:, 0] -= 1010922
xyz[:, 1] -= 6810940


# x = np.linspace(xmin-2*Lh, xmax+2*Lh)
# y = np.linspace(ymin-2*Lh, ymax+2*Lh)

# centers = np.meshgrid(x, y)
# z = 500 + 1000 * np.sin((centers[0]-xmin) * (2*np.pi)/(xmax-xmin))

# xyz = np.transpose(np.vstack([a.ravel() for a in centers] +
#                             [z.ravel(),]))

dtm = urg.elevation_surface(xyz)
MT.to_vtu(MT.TSurf.make(*dtm.as_arrays()), "topography.vtu")

# print('isotropic remeshing')
# mesh.isotropic_remeshing(3000.)
# mesh.add_surface(dtm)

# dump_model(mesh, 'model_with_dtm.vtu')

# from matplotlib import pyplot as plt

# plt.clf()
# plt.imshow(data)
# plt.show()
