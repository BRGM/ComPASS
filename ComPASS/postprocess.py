# -*- coding: utf-8 -*-
#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import sys
import shutil
import glob, os, re
from pathlib import Path
import click
import numpy as np
from .dumps import Dumper
from .dump_wells import _well_vtu_filename as well_vtu
from .utils import create_directories
from .utils import units
import vtkwriters as vtkw
from .utils.adjacencies import filter_adjacency_table


class MeshDistribution:
    def __init__(self, filename):
        def extract_integer_list(s):
            match = re.match("\[(.*)\]", s.strip())
            assert match is not None
            return [int(si.strip()) for si in match.group(1).split()]

        vars = {
            "nb_procs": lambda s: int(s),
            "nb_own_cells": extract_integer_list,
            "nb_own_nodes": extract_integer_list,
            "nb_own_faces": extract_integer_list,
            "nb_own_fractures": extract_integer_list,
        }
        for name in vars:
            setattr(self, name, None)
        with open(filename) as f:
            for line in f:
                l = line.lower().strip()
                for name, extract in vars.items():
                    if getattr(self, name) is None:
                        match = re.match(name + "\s*=\s*(.*)", l)
                        if match is not None:
                            setattr(self, name, extract(match.group(1)))
                            break

    def __str__(self):
        return "\n".join(
            [
                "%s = %d" % (s, getattr(self, s))
                for s in (
                    "nb_procs",
                    "nb_own_cells",
                    "nb_own_nodes",
                    "nb_own_faces",
                    "nb_own_fractures",
                )
            ]
        )


class PostProcessor:
    def __init__(self, directory, time_unit):
        # No simulation object during postprocess -> None
        self.dumper = Dumper(None, directory)
        self.paraview_directory = os.path.join(
            self.dumper.to_output_directory(), "paraview"
        )
        create_directories(self.paraview_directory)
        self.vtu_directory = os.path.join(self.paraview_directory, "vtu")
        create_directories(self.vtu_directory)
        self.pvwells_directory = (
            self.vtu_directory
        )  # everything in the same vtu directory
        create_directories(self.pvwells_directory)
        self.distribution = MeshDistribution(
            self.dumper.to_output_directory("own_elements")
        )
        self.vertices_type = None
        self.data_types = None
        self.proc_id_type = np.dtype("i")
        assert time_unit > 0
        self.time_unit = time_unit

    def to_paraview_directory(self, filename=""):
        return os.path.join(self.paraview_directory, filename)

    def to_pvwells_directory(self, filename=""):
        return os.path.join(self.pvwells_directory, filename)

    def to_vtu_directory(self, filename=""):
        return os.path.join(self.vtu_directory, filename)

    def mesh_filename(self, proc):
        assert proc < self.distribution.nb_procs
        return self.dumper.mesh_filename(proc)

    def states_filename(self, proc, tag=""):
        assert proc < self.distribution.nb_procs
        return self.dumper.states_filename(proc, tag)

    def write_mesh_vtu(
        self,
        proc,
        basename,
        own_only=True,
        dump_procid=False,
        nodedata=None,
        celldata=None,
        fracdata=None,
    ):
        assert proc < self.distribution.nb_procs
        nb_own_cells = self.distribution.nb_own_cells[proc]
        nb_own_fractures = self.distribution.nb_own_fractures[proc]
        if dump_procid:
            own_only = True
            assert nodedata is None and celldata is None and fracdata is None
            celldata = {
                "proc": np.array(np.tile(proc, nb_own_cells), dtype=self.proc_id_type)
            }
            if nb_own_fractures > 0:
                fracdata = {
                    "proc": np.array(
                        np.tile(proc, nb_own_fractures), dtype=self.proc_id_type
                    )
                }
        if nodedata is None:
            nodedata = {}
        if celldata is None:
            celldata = {}
        if fracdata is None:
            fracdata = {}
        meshfile = self.mesh_filename(proc)
        assert os.path.exists(meshfile)
        mesh = np.load(meshfile)
        if self.vertices_type is None:
            self.vertices_type = mesh["vertices"].dtype
        assert self.vertices_type == mesh["vertices"].dtype
        proc_label = self.dumper.proc_label(proc)
        piecefile = self.to_vtu_directory("%s_%s.vtu" % (basename, proc_label))
        coc_to_list = lambda s: np.split(mesh[f"{s}_values"], mesh[f"{s}_offsets"][:-1])
        cellfaces = coc_to_list("cellfaces")
        facenodes = coc_to_list("facenodes")
        cells = [[facenodes[face] for face in faces] for faces in cellfaces]
        if own_only:
            cells = cells[:nb_own_cells]
        vtkw.write_vtu(
            vtkw.polyhedra_vtu_doc(
                mesh["vertices"],
                cells,
                pointdata=nodedata,
                celldata=celldata,
            ),
            piecefile,
        )
        # FIXME: we could optimize selecting only fracture nodes and skipping empty files
        if fracdata:
            # will triger a TypeError if array sizes are not the same
            fracdata_size = np.unique([a.shape[0] for a in fracdata.values()])
            assert len(fracdata_size) == 1
            fracdata_size = int(fracdata_size[0])
            if fracdata_size > 0:
                fracpiecefile = self.to_vtu_directory(
                    "fracture_%s_%s.vtp" % (basename, proc_label)
                )
                cell_nodes_offsets = mesh["fracturenodes_offsets"]
                cell_nodes = mesh["fracturenodes_values"]
                cell_types = mesh["fracture_types"]
                if own_only:
                    if nb_own_fractures == 0:
                        return piecefile, None
                    cell_nodes_offsets = cell_nodes_offsets[:nb_own_fractures]
                    cell_nodes = cell_nodes[: cell_nodes_offsets[-1]]
                # FIXME: this is suboptimal because we should use the underlying structures
                faces_nodes = np.split(cell_nodes, cell_nodes_offsets[:-1])
                used, faces_nodes = filter_adjacency_table(faces_nodes)
                vtkw.write_vtp(
                    vtkw.vtp_doc(
                        mesh["vertices"][used],
                        faces_nodes,
                        celldata=fracdata,
                    ),
                    fracpiecefile,
                )
                return piecefile, fracpiecefile
        return piecefile, None

    def collect_proc_ids(self):
        pieces, fracpieces = [], []
        for proc in range(self.distribution.nb_procs):
            piecefiles = self.write_mesh_vtu(proc, "mesh", dump_procid=True)
            assert piecefiles[0] is not None
            pieces.append(os.path.relpath(piecefiles[0], self.to_paraview_directory()))
            if piecefiles[1] is not None:
                fracpieces.append(
                    os.path.relpath(piecefiles[1], self.to_paraview_directory())
                )
        vtkw.write_pvtu(
            vtkw.pvtu_doc(
                self.vertices_type,
                pieces,
                celldata_types={"proc": self.proc_id_type},
            ),
            self.to_paraview_directory("mesh.pvtu"),
        )
        if len(fracpieces) > 0:
            vtkw.write_pvtp(
                vtkw.pvtp_doc(
                    self.vertices_type,
                    fracpieces,
                    celldata_types={"proc": self.proc_id_type},
                ),
                self.to_paraview_directory("fractures_mesh"),
            )

    def collect_snapshots(self):
        snapshots_file = self.dumper.to_output_directory("snapshots")
        assert os.path.exists(snapshots_file)
        with open(snapshots_file) as f:
            snapshots = {}
            for line in f:

                l = line.strip().split()
                if len(l) >= 2:
                    snapshots[l[0]] = float(l[1])
        return snapshots

    def extract_data(self, proc, tag):
        datafile = self.states_filename(proc, tag)
        assert os.path.exists(datafile)
        data = np.load(datafile)

        def extract_at(location):
            return {
                name.replace(location + "_", ""): data[name]
                for name in data.files
                if name.startswith(location)
            }

        data_arrays = {
            location: extract_at(location) for location in ("node", "cell", "fracture")
        }
        datatype_checks = {
            location: {
                name: property.dtype
                if len(property.shape) == 1
                else (property.dtype, property.shape[1])
                for name, property in properties.items()
            }
            for location, properties in data_arrays.items()
        }
        if self.data_types is None:
            self.data_types = datatype_checks
        assert self.data_types == datatype_checks
        return data_arrays

    def collect_states(self, convert_temperature):
        snapshots = self.collect_snapshots()
        pvd = {}
        for tag in snapshots:
            all_pieces = []
            for proc in range(self.distribution.nb_procs):
                data_arrays = self.extract_data(proc, tag)
                if convert_temperature:
                    for location in data_arrays.keys():
                        for property in data_arrays[location].keys():
                            if "temperature" in property:
                                data_arrays[location][property] -= 273.15
                pieces = self.write_mesh_vtu(
                    proc,
                    "states_" + tag,
                    nodedata=data_arrays["node"],
                    celldata=data_arrays["cell"],
                    fracdata=data_arrays["fracture"],
                )
                all_pieces.append(pieces)
            pvtufile = self.to_vtu_directory("state_" + tag + ".pvtu")
            vtkw.write_pvtu(
                vtkw.pvtu_doc(
                    self.vertices_type,
                    [
                        os.path.relpath(pieces[0], self.to_vtu_directory())
                        for pieces in all_pieces
                    ],
                    pointdata_types=self.data_types["node"],
                    celldata_types=self.data_types["cell"],
                ),
                pvtufile,
            )
            pvtus = (pvtufile,)
            if any(pieces[1] is not None for pieces in all_pieces):
                fracpvtpfile = self.to_vtu_directory("fracture_state_" + tag + ".pvtp")
                vtkw.write_pvtp(
                    vtkw.pvtp_doc(
                        self.vertices_type,
                        [
                            os.path.relpath(pieces[1], self.to_vtu_directory())
                            for pieces in all_pieces
                            if pieces[1] is not None
                        ],
                        celldata_types=self.data_types["fracture"],
                    ),
                    fracpvtpfile,
                )
                pvtus = (pvtufile, fracpvtpfile)
            pvd[snapshots[tag] / self.time_unit] = pvtus
        vtkw.write_pvd(
            vtkw.pvd_doc(
                [
                    (t, os.path.relpath(pvtus[0], self.to_paraview_directory()))
                    for t, pvtus in pvd.items()
                ]
            ),
            self.to_paraview_directory("states"),
        )
        if any(len(pvtus) > 1 for pvtus in pvd.values()):
            assert all(len(pvtus) > 1 for pvtus in pvd.values())
            vtkw.write_pvd(
                vtkw.pvd_doc(
                    [
                        (t, os.path.relpath(pvtus[1], self.to_paraview_directory()))
                        for t, pvtus in pvd.items()
                    ]
                ),
                self.to_paraview_directory("fracture_states"),
            )

    def collect_wells(self, convert_temperature):
        # FIXME: temperature conversion is not available (cf. gitlab issue #205)
        well_ids_file = Path(self.dumper.to_output_directory("well_ids"))
        if not well_ids_file.exists():
            return
        with open(well_ids_file) as f:
            if len("".join(f.readlines()).strip()) == 0:
                return
        well_ids = np.loadtxt(well_ids_file, ndmin=1, dtype=np.uint64)
        snapshots = self.collect_snapshots()
        for well in well_ids:
            pvd = {}
            for tag in snapshots:
                vtu = self.dumper.to_wells_directory(well_vtu(well, tag))
                assert os.path.exists(vtu)
                # copy file in paraview directory
                new_vtu = self.to_pvwells_directory(os.path.basename(vtu))
                shutil.copyfile(vtu, new_vtu)
                pvd[snapshots[tag] / self.time_unit] = new_vtu
            vtkw.write_pvd(
                vtkw.pvd_doc(
                    [
                        (t, os.path.relpath(vtu, self.to_paraview_directory()))
                        for t, vtu in pvd.items()
                    ]
                ),
                self.to_paraview_directory(f"well-{well:04d}"),
            )


def postprocess(
    directory,
    collect_procs_id=False,
    collect_states=True,
    convert_temperature=True,
    collect_wells=True,
    time_unit="year",
):
    """postprocess a set of directories where output from ComPASS simulations are stored (typically something like output-scriptname)

    :param directory: directory to process
    :param collect_procs_id: boolean flag to collect procs ids and output mesh partitioning
    :param collect_states: boolean flag to collect physical states
    :param convert_temperature: boolean flag to convert Kelvin to Celsius degrees
    :param collect_wells: boolean flag to collect well information
    :param time_unit: a string among second, minute, hour, day, year, defaults to year

    """
    directory = Path(directory)
    print("processing results in", directory)
    something_done = False
    time_unit_value = {
        "second": 1,
        "minute": units.minute,
        "hour": units.hour,
        "day": units.day,
        "year": units.year,
    }[time_unit]
    if directory.is_dir() and (collect_procs_id or collect_states or collect_wells):
        pp = PostProcessor(directory, time_unit_value)
        if collect_procs_id:
            pp.collect_proc_ids()
            something_done = True
        if collect_states:
            pp.collect_states(convert_temperature)
            something_done = True
        if collect_wells:
            if not convert_temperature:
                print("***************** WARNING ******************")
                print()
                print("Well temperature is output in Celsius degree")
                print("          (cf. gitlab issue #205)")
                print()
                print("********************************************")
            pp.collect_wells(convert_temperature)
            something_done = True
    if not something_done:
        print()
        print("************* WARNING **************")
        print()
        print(f"Nothing done in: {directory}")
        print()
        print("************************************")
        print()


@click.command(name="postprocess")
@click.option(
    "-p",
    "--procs",
    "collect_procs_id",
    is_flag=True,
    default=False,
    help="ouput paraview/mesh.pvtu file with cell distribution (proc variable)",
)
@click.option(
    "-s",
    "--states",
    "collect_states",
    is_flag=True,
    default=False,
    help="ouput paraview files with transient states as pvd files to be found in the path/paraview directory",
)
@click.option(
    "-w",
    "--wells",
    "collect_wells",
    is_flag=True,
    default=False,
    help="ouput wells transient states in paraview pvd files to be found in the path/paraview directory",
)
@click.option(
    "-C",
    "--Celsius",
    "convert_temperature",
    is_flag=True,
    default=False,
    help="convert temperature from Kelvin to Celsius degrees",
)
@click.option(
    "-t",
    "--time-unit",
    "time_unit",
    default="year",
    type=click.Choice(
        ["second", "minute", "hour", "day", "year"], case_sensitive=False
    ),
    show_default=True,
    help="select a time unit",
)
@click.argument("directories", nargs=-1)
def postprocess_command(
    collect_procs_id,
    collect_states,
    collect_wells,
    convert_temperature,
    time_unit,
    directories,
):
    """postprocess a set of directories where output from ComPASS simulations are stored (typically something like output-scriptname)"""
    for directory in directories:
        postprocess(
            directory,
            collect_procs_id,
            collect_states,
            convert_temperature,
            collect_wells,
            time_unit,
        )


if __name__ == "__main__":
    postprocess_command(prog_name="compass postprocess")
