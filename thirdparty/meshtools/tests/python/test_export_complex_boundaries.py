import numpy as np
import MeshTools as MT
from MeshTools.CGALWrappers import mesh_implicit_domains_boundaries

as_tsurf = lambda t: MT.TSurf.make(t[0], MT.idarray(t[1]))

def export_meshes(meshes, basename="mesh"):
    for i, mesh in enumerate(meshes):
        MT.to_vtu(
            as_tsurf(mesh.as_arrays()),
            '%s%02d.vtu' % (basename, i)
        )

export_meshes(mesh_implicit_domains_boundaries(),
              "implicit_domains_boundaries-")
