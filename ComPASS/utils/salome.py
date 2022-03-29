# -*- coding: utf-8 -*-

from pathlib import Path
from collections import defaultdict, OrderedDict, namedtuple
import importlib
import numpy as np
import MeshTools as MT
import vtkwriters as vtkw
import ComPASS.mpi as mpi

# TODO: dict of mask values could be broadcasted using mpi4py

_locations = ["nodes", "faces", "cells"]
ElementValues = namedtuple(
    "ElementValues",
    _locations,
    defaults=[
        None,
    ]
    * len(_locations),
)


class IdProvider:
    def __init__(self, nodes_provider):
        self._nodes_provider = nodes_provider
        self._id_map = None

    @property
    def id_map(self):
        if self._id_map is None:
            eltnodes = np.copy(self._nodes_provider())
            eltnodes.sort(axis=-1)
            self._id_map = {tuple(nodes): k for k, nodes in enumerate(eltnodes)}
        return self._id_map

    def __getitem__(self, nodes):
        return self.id_map[tuple(sorted(nodes))]


class MTWrapper:
    def __init__(self, mesh):
        super().__init__()
        self._mesh = mesh
        self._cellnodes = None
        self._cell_id = IdProvider(lambda: self.cellnodes)

    @property
    def vertices(self):
        return self._mesh.vertices_array()

    @property
    def facenodes(self):
        return self._mesh.connectivity.faces.nodes

    @property
    def cellnodes(self):
        if self._cellnodes is None:
            offsets, cellnodes = self._mesh.cells_nodes_as_COC()
            cellnodes = cellnodes.array_view()
            cellnodes.shape = -1, 4
            self._cellnodes = cellnodes
        return self._cellnodes

    @property
    def nb_vertices(self):
        return self._mesh.nb_vertices

    @property
    def nb_faces(self):
        return self._mesh.nb_faces

    @property
    def nb_cells(self):
        return self._mesh.nb_cells

    def face_id(self, nodes):
        return self._mesh.face_id(MT.Triangle(nodes))

    def cell_id(self, nodes):
        return self._cell_id[nodes]

    def to_vtu(self, filepath, pointdata={}, celldata={}):
        MT.to_vtu(self._mesh, str(filepath), pointdata=pointdata, celldata=celldata)


class CompassWrapper:
    def __init__(self, simulation):
        self._simulation = simulation
        self._facenodes = None
        self._cellnodes = None
        self._face_id = IdProvider(lambda: self.facenodes)
        self._cell_id = IdProvider(lambda: self.cellnodes)

    @property
    def vertices(self):
        return self._simulation.vertices()

    @property
    def facenodes(self):
        if self._facenodes is None:
            connectivity = self._simulation.get_connectivity().NodebyFace
            facenodes = np.copy(connectivity.contiguous_content())
            facenodes -= 1  # Fortran indexing -> C indexing
            offsets = connectivity.offsets()
            assert np.all(
                offsets[1:] - offsets[:-1] == 3
            ), "All facets should be triangles!"
            facenodes.shape = -1, 3
            self._facenodes = facenodes
        return self._facenodes

    @property
    def cellnodes(self):
        if self._cellnodes is None:
            connectivity = self._simulation.get_connectivity().NodebyCell
            cellnodes = np.copy(connectivity.contiguous_content())
            cellnodes -= 1  # Fortran indexing -> C indexing
            offsets = connectivity.offsets()
            assert np.all(
                offsets[1:] - offsets[:-1] == 4
            ), "All cells should be tetraedra!"
            cellnodes.shape = -1, 4
            self._cellnodes = cellnodes
        return self._cellnodes

    @property
    def nb_vertices(self):
        return self.vertices.shape[0]

    @property
    def nb_faces(self):
        connectivity = self._simulation.get_connectivity().NodebyFace
        return connectivity.offsets().shape[0] - 1

    @property
    def nb_cells(self):
        connectivity = self._simulation.get_connectivity().NodebyCell
        return connectivity.offsets().shape[0] - 1

    def face_id(self, nodes):
        return self._face_id[nodes]

    def cell_id(self, nodes):
        return self._cell_id[nodes]

    def to_vtu(self, filepath, pointdata={}, celldata={}):
        vtkw.write_vtu(
            vtkw.vtu_doc(
                self.vertices, self.cellnodes, pointdata=pointdata, celldata=celldata
            ),
            str(filepath),
        )


class SalomeGroupData:
    def __init__(self, ids=None):
        self.ids = ids
        self.mask = None


class SalomeGroup:
    _cells_id_map = defaultdict(lambda: None)

    def __init__(
        self, mesh, nodes=None, faces=None, cells=None, facenodes=None, cellnodes=None
    ):
        self._mesh = mesh
        self._nodes = SalomeGroupData(nodes)
        self._faces = SalomeGroupData()
        self._cells = SalomeGroupData()
        if facenodes is not None:
            assert (
                nodes is None and cells is None and cellnodes is None and faces is None
            )
            self._add_facenodes(facenodes)
        if cellnodes is not None:
            assert (
                nodes is None and faces is None and facenodes is None and cells is None
            )
            self._add_cellnodes(cellnodes)
        if faces is not None:
            assert (
                nodes is None
                and cells is None
                and cellnodes is None
                and facenodes is None
            )
            self._add_faces(faces)
        if cells is not None:
            assert (
                nodes is None
                and faces is None
                and facenodes is None
                and cellnodes is None
            )
            self._add_cells(cells)

    @property
    def has_no_ids(self):
        return (
            self._nodes.ids is None
            and self._faces.ids is None
            and self._cells.ids is None
        )

    def _add_facenodes(self, faces):
        assert self.has_no_ids
        self._nodes.ids = np.unique(faces)
        mesh = self._mesh
        self._faces.ids = np.array([mesh.face_id(nodes) for nodes in faces])
        assert np.all(
            self._faces.ids < mesh.nb_faces
        ), f"A face from {name} group could not be found."

    def _add_faces(self, faces):
        assert self.has_no_ids
        self._faces.ids = faces
        mesh = self._mesh
        if mesh is not None:
            facenodes = mesh.facenodes
            self._nodes.ids = np.unique(
                np.array([np.array(facenodes[fk]) for fk in faces])
            )

    def _add_cellnodes(self, cells):
        assert self.has_no_ids
        self._nodes.ids = np.unique(cells)
        mesh = self._mesh
        self._cells.ids = np.array([mesh.cell_id(nodes) for nodes in cells])

    def _add_cells(self, cells):
        assert self.has_no_ids
        self._cells.ids = cells
        mesh = self._mesh
        if mesh is not None:
            cellnodes = mesh.cellnodes
            self._nodes.ids = np.unique(
                np.array([np.array(cellnodes[ck]) for ck in cells])
            )

    @property
    def nodes(self):
        # No conversion implemtented yet
        return self._nodes.ids

    @property
    def faces(self):
        # No conversion implemtented yet
        return self._faces.ids

    @property
    def cells(self):
        # No conversion implemtented yet
        return self._cells.ids


def _extract_groups_module_info(groups_module="GROUPS"):
    groups = importlib.import_module(groups_module)
    group_names = [name for name in dir(groups) if not name.startswith("__")]
    return {name: MT.idarray(getattr(groups, name)) - 1 for name in group_names}


class SalomeMeshInfo:
    def __init__(self, mesh=None):
        self._dict = OrderedDict()
        self._mesh = mesh

    def _build_from_groups(self, groups):
        mesh = self._mesh
        for name, vids in groups.items():
            if vids.ndim == 1:  # nodes
                self._dict[name] = SalomeGroup(mesh, nodes=vids)
            elif vids.ndim == 2 and vids.shape[1] == 3:  # triangles
                self._dict[name] = SalomeGroup(mesh, facenodes=vids)
            elif vids.ndim == 2 and vids.shape[1] == 4:  # tets
                self._dict[name] = SalomeGroup(mesh, cellnodes=vids)
            else:
                assert False, f"Could not handle group: {name}"

    def build_from_groups_module(self, groups_module="GROUPS"):
        groups = _extract_groups_module_info(groups_module)
        self._build_from_groups(groups)

    def build_from_npz(self, filename):
        groups = np.load(filename)
        self._build_from_groups(groups)

    def __getattr__(self, name):
        try:
            return self.__dict__[name]
        except:
            pass
        try:
            return self._dict[name]
        except KeyError:
            return AttributeError(f"No group with name {name} found!")

    def groups(self):
        return self._dict.keys()

    def get(self, key):
        return self._dict.get(key)

    def items(self):
        return self._dict.items()

    @property
    def mesh(self):
        return self._mesh

    def to_vtu_block(self, filename):
        filepath = Path(filename)
        filepath.parent.mkdir(parents=True, exist_ok=True)
        mesh = self._mesh
        assert mesh is not None
        pointdata = {}
        celldata = {}
        for name in self.groups():
            group = getattr(self, name)
            if group.nodes is not None:
                where = np.zeros(mesh.nb_vertices, dtype=np.int32)
                where[group.nodes] = 1
                pointdata[name] = where
            if group.cells is not None:
                where = np.zeros(mesh.nb_cells, dtype=np.int32)
                where[group.cells] = 1
                celldata[name] = where
        mesh.to_vtu(filepath, pointdata=pointdata, celldata=celldata)

    def faces_to_multiblock(self, filename):
        filepath = Path(filename)
        filepath.parent.mkdir(parents=True, exist_ok=True)
        subdir = filepath.parent / f"{filepath.stem}-elements"
        subdir.mkdir(exist_ok=True)
        mesh = self._mesh
        vertices = mesh.vertices
        assert mesh is not None
        elements = []
        for name in self.groups():
            group = getattr(self, name)
            if group.faces is not None:
                facenodes = mesh.facenodes
                group_faces = np.array([np.array(facenodes[fk]) for fk in group.faces])
                group_vertices, group_faces = np.unique(
                    group_faces, return_inverse=True
                )
                assert np.all(group_vertices == group.nodes)
                tsurf = MT.TSurf.make(vertices[group_vertices], MT.idarray(group_faces))
                element_path = subdir / f"{name}.vtu"
                MT.to_vtu(tsurf, str(element_path))
                elements.append((name, str(element_path)))
        vtkw.write_vtm(vtkw.vtm_doc(elements), str(filepath))

    def _compute_node_masks(self):
        mask = 1
        for name, group in self.items():
            if group.nodes is not None:
                if group.faces is None and group.cells is None:
                    group._nodes.mask = mask
                    mask *= 2
                else:
                    group._nodes.mask = None

    def _compute_face_masks(self):
        mask = 1
        for name, group in self.items():
            if group.faces is not None:
                if group.cells is None:
                    group._faces.mask = mask
                    mask *= 2
                else:
                    group._faces.mask = None

    def _compute_cell_masks(self):
        mask = 1
        for name, group in self.items():
            if group.cells is not None:
                group._cells.mask = mask
                mask *= 2
            else:
                group._cells.mask = None

    def compute_masks(self):
        self._compute_node_masks()
        self._compute_face_masks()
        self._compute_cell_masks()

    def _add_node_mask(self, masks, new_name, new_group):
        # find highest mask in nodes group
        mask = 1
        for name, group in masks.items():
            if group.nodes is not None:
                mask = max(mask, 2 * group.nodes)
        masks[new_name] = ElementValues(mask, None, None)

    def _add_face_mask(self, masks, new_name, new_group):
        # find highest mask in faces group
        mask = 1
        for name, group in masks.items():
            if group.faces is not None:
                mask = max(mask, 2 * group.faces)
        masks[new_name] = ElementValues(None, mask, None)

    def _add_cell_mask(self, masks, new_name, new_group):
        # find highest mask in cells group
        mask = 1
        for name, group in masks.items():
            if group.cells is not None:
                mask = max(mask, 2 * group.cells)
        masks[new_name] = ElementValues(None, None, mask)

    def add_mask(self, masks, new_name, new_group):
        # check if new_name does not already exist
        assert (
            masks.get(new_name) is None
        ), "Cannot give twice the same name in groups: {new_name}"
        # creates a SalomeGroup from new_group
        mesh = self._mesh
        if new_group.ndim == 1:  # nodes
            sgroup = SalomeGroup(mesh, nodes=new_group)
            self._dict[new_name] = sgroup
            self._add_node_mask(masks, new_name, sgroup)
        elif new_group.ndim == 2 and new_group.shape[1] == 3:  # triangles
            sgroup = SalomeGroup(mesh, facenodes=new_group)
            self._dict[new_name] = sgroup
            self._add_face_mask(masks, new_name, sgroup)
        elif new_group.ndim == 2 and new_group.shape[1] == 4:  # tets
            sgroup = SalomeGroup(mesh, cellnodes=new_group)
            self._dict[new_name] = sgroup
            self._add_cell_mask(masks, new_name, sgroup)
        else:
            assert False, f"Could not handle group: {new_name}"

    def collect_masks(self):
        return OrderedDict(
            [
                (
                    name,
                    ElementValues(
                        group._nodes.mask, group._faces.mask, group._cells.mask
                    ),
                )
                for name, group in self.items()
            ]
        )

    def compute_and_collect_masks(self):
        self.compute_masks()
        return self.collect_masks()

    def add_and_collect_masks(self, masks, name, group):
        self.add_mask(masks, name, group)
        return masks

    def set_flags(self, flags, recompute_masks=True):
        if recompute_masks:
            self.compute_masks()
        for name, group in self.items():
            for location in _locations:
                where = getattr(group, location)
                if where is not None:
                    mask = getattr(group, f"_{location}").mask
                    if mask is not None:
                        location_flags = getattr(flags, location)
                        maxmask = 2 ** (8 * location_flags.dtype.itemsize)
                        assert (
                            mask <= maxmask
                        ), "Cannot encode more than {maxmask} groups!"
                        location_flags[where] |= mask

    def rebuild_from_flags(self, flags, masks, mesh=None):
        self._mesh = mesh
        self._dict = OrderedDict()
        for name, values in masks.items():
            assert list(values).count(None) == 2
            for mask, location in zip(values, _locations):
                if mask is not None:
                    where = np.nonzero(getattr(flags, location) & mask)[0]
                    if where.size == 0:
                        self._dict[name] = SalomeGroup(mesh)
                    else:
                        self._dict[name] = SalomeGroup(mesh, **{location: where})


def _load_mesh_elements(nodes_file, tets_file):
    vertices = np.loadtxt(nodes_file, usecols=(1, 2, 3))
    tets = np.loadtxt(tets_file, dtype=np.ulonglong, usecols=(1, 2, 3, 4))
    assert np.all(tets > 0)
    cells = MT.idarray(tets - 1)  # Salome indexing starts at 1
    return vertices, cells


def _build_mesh(nodes_file, tets_file):
    return MT.TetMesh.make(*_load_mesh_elements(nodes_file, tets_file))


def load(nodes_file="NODES.txt", tets_file="TETRAS.txt", groups_module="GROUPS"):
    mesh = _build_mesh(nodes_file, tets_file)
    mesh_info = SalomeMeshInfo(MTWrapper(mesh))
    mesh_info.build_from_groups_module(groups_module)
    return mesh, mesh_info


def _groups_npz(basename):
    return basename + ".groups"


def compress(
    basename,
    nodes_file="NODES.txt",
    tets_file="TETRAS.txt",
    groups_module="GROUPS",
    verbose=False,
):
    if verbose:
        print("Extracting mesh elements")
    vertices, cells = _load_mesh_elements(nodes_file, tets_file)
    if verbose:
        print(f"{vertices.shape[0]} vertices and {cells.shape[0]} cells were extracted")
    np.savez(basename, vertices=vertices, cells=cells)
    if verbose:
        print("Extracting groups")
    groups = _extract_groups_module_info(groups_module)
    if verbose:
        print("Extracted groups:", ", ".join(groups.keys()))
    np.savez(_groups_npz(basename), **groups)


def load_compressed(basename):
    data = np.load(basename + ".npz")
    mesh = MT.TetMesh.make(data["vertices"], data["cells"])
    mesh_info = SalomeMeshInfo(MTWrapper(mesh))
    mesh_info.build_from_npz(_groups_npz(basename) + ".npz")
    return mesh, mesh_info


class SalomeWrapper:
    def __init__(
        self,
        simulation,
        nodes_file="NODES.txt",
        tets_file="TETRAS.txt",
        groups_module="GROUPS",
        verbose=True,
        compressed_basename=None,
    ):
        self._simulation = simulation
        if mpi.is_on_master_proc:
            if compressed_basename is None:
                mesh, info = load(nodes_file, tets_file, groups_module)
            else:
                mesh, info = load_compressed(compressed_basename)
            masks = info.compute_and_collect_masks()
            if verbose:
                print("Available Salome groups:", ", ".join(info.groups()))
            self.mesh = mesh
            self.info = info
        else:
            self.mesh = None
            self.info = SalomeMeshInfo()
            masks = None
        self.masks = mpi.communicator().bcast(masks, root=mpi.master_proc_rank)

    @property
    def flags_setter(self, recompute_masks=True):
        simulation = self._simulation

        def setter():
            self.info.set_flags(
                ElementValues(
                    simulation.global_nodeflags(),
                    simulation.global_faceflags(),
                    simulation.global_cellflags(),
                ),
                recompute_masks,
            )

        return setter

    def rebuild_locally(self):
        self.mesh = CompassWrapper(self._simulation)
        simulation = self._simulation
        flags = ElementValues(
            simulation.nodeflags(),
            simulation.faceflags(),
            simulation.cellflags(),
        )
        self.info.rebuild_from_flags(flags, self.masks, self.mesh)

    def add_group(self, name, group, verbose=True):
        if mpi.is_on_master_proc:
            masks = self.info.add_and_collect_masks(self.masks, name, group)
            if verbose:
                print("Available Salome groups:", ", ".join(self.info.groups()))
        else:
            masks = None
        self.masks = mpi.communicator().bcast(masks, root=mpi.master_proc_rank)
