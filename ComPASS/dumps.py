#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import os
import numpy as np

from . import mpi
from . import dump_wells as dw
from .utils import create_directories
from .runtime import to_output_directory


class Dumper:
    def __init__(self, simulation, output_directory=None):
        self.simulation = simulation
        if output_directory is None:
            output_directory = to_output_directory("")
        self.output_directory = os.path.abspath(output_directory)
        self.mesh_directory = os.path.join(self.output_directory, "mesh")
        self.states_directory = os.path.join(self.output_directory, "states")
        self.wells_directory = os.path.join(self.output_directory, "wells")
        self.simulation_running = False
        self.proc_label_pattern = (
            "proc_{:06d}"  # IMPROVE: hard coded maximum nb of procs
        )
        create_directories(self.to_output_directory())
        create_directories(self.to_mesh_directory())
        create_directories(self.to_states_directory())
        create_directories(self.to_wells_directory())

    def to_output_directory(self, filename=""):
        return os.path.join(self.output_directory, filename)

    def to_mesh_directory(self, filename=""):
        return os.path.join(self.mesh_directory, filename)

    def to_states_directory(self, filename=""):
        return os.path.join(self.states_directory, filename)

    def to_wells_directory(self, filename=""):
        return os.path.join(self.wells_directory, filename)

    def proc_label(self, proc):
        return self.proc_label_pattern.format(proc)

    def mesh_filename(self, proc):
        return self.to_mesh_directory("mesh_%s.npz" % self.proc_label(proc))

    def states_filename(self, proc, tag=""):
        if tag:
            tag += "_"
        return self.to_states_directory("state_%s%s.npz" % (tag, self.proc_label(proc)))

    def start_simulation(self):
        assert not self.simulation_running
        self.simulation_running = True
        self.dump_own_element_numbers()
        self.dump_wellids()
        self.dump_mesh()

    @mpi.on_master_proc
    def dump_own_element_numbers(self):
        filename = self.to_output_directory("own_elements")
        communicator = mpi.communicator()
        np_linewidth_backup = np.get_printoptions()["linewidth"]
        np.set_printoptions(
            linewidth=np.inf
        )  # we want all outputs below on a single line
        with open(filename, "w") as f:
            print("# Number of procs", file=f)
            print("nb_procs =", communicator.size, file=f)
            print("nb_own_cells =", self.simulation.nb_cells_own(), file=f)
            print("nb_own_nodes =", self.simulation.nb_nodes_own(), file=f)
            print("nb_own_faces =", self.simulation.nb_faces_own(), file=f)
            print("nb_own_fractures =", self.simulation.nb_fractures_own(), file=f)
        np.set_printoptions(linewidth=np_linewidth_backup)

    def dump_wellids(self):
        well_ids = dw.collect_well_ids(self.simulation)
        well_ids = mpi.communicator().gather(well_ids, root=mpi.master_proc_rank)
        if mpi.is_on_master_proc:
            filename = self.to_output_directory("well_ids")
            with open(filename, "w") as f:
                for ids in well_ids:
                    for i in ids:
                        print(i, file=f)

    def dump_mesh(self):
        connectivity = self.simulation.get_connectivity()
        fracture_nodes_coc = self.simulation.get_nodes_by_fractures()
        fracture_nodes = [
            np.array(nodes) - 1 for nodes in fracture_nodes_coc
        ]  # switch first node indexing from 1 to 0
        fracturenodes_offsets = np.cumsum([len(a) for a in fracture_nodes])
        fracturenodes_values = (
            np.hstack(fracture_nodes) if len(fracture_nodes) > 0 else np.array([])
        )
        fracture_faces = self.simulation.frac_face_id()
        fracture_types = self.simulation.facetypes()[
            fracture_faces - 1
        ]  # switch first node indexing from 1 to 0
        np.savez(
            self.mesh_filename(mpi.proc_rank),
            vertices=self.simulation.vertices().view(dtype=np.double).reshape((-1, 3)),
            cellnodes_offsets=connectivity.NodebyCell.offsets()[
                1:
            ],  # VTK does not use the first 0 offset
            cellnodes_values=connectivity.NodebyCell.contiguous_content()
            - 1,  # switch first node indexing from 1 to 0
            celltypes=self.simulation.celltypes(),
            fracturenodes_offsets=fracturenodes_offsets,
            fracturenodes_values=fracturenodes_values,
            fracture_types=fracture_types,
        )

    def dump_states(self, tag="", dump_fluxes=True):
        assert self.simulation_running
        node_states = self.simulation.node_states()
        cell_states = self.simulation.cell_states()
        fracture_states = self.simulation.fracture_states()
        dumped_states = {
            "node_pressure": node_states.p,
            "node_temperature": node_states.T,
            "cell_pressure": cell_states.p,
            "cell_temperature": cell_states.T,
            "fracture_pressure": fracture_states.p,
            "fracture_temperature": fracture_states.T,
        }
        nbphases = self.simulation.number_of_phases()
        phase_names = [
            name
            for name in self.simulation.Phase.__dict__.keys()
            if not name.startswith("__")
        ]
        nbcomponents = self.simulation.number_of_components()
        comp_names = [
            name
            for name in self.simulation.Component.__dict__.keys()
            if not name.startswith("__")
        ]
        for phase in range(nbphases):
            phasename = phase_names[phase]
            dumped_states["cell_saturation_phase_%s" % (phasename)] = cell_states.S[
                :, phase
            ]
            dumped_states["node_saturation_phase_%s" % (phasename)] = node_states.S[
                :, phase
            ]
            dumped_states[
                "fracture_saturation_phase_%s" % (phasename)
            ] = fracture_states.S[:, phase]
        for phase in range(nbphases):
            phasename = phase_names[phase]
            if nbcomponents > 1:
                for comp in range(nbcomponents):
                    compname = comp_names[comp]
                    dumped_states[
                        "cell_comp_%s_in_phase_%s" % (compname, phasename)
                    ] = cell_states.C[:, phase, comp]
                    dumped_states[
                        "node_comp_%s_in_phase_%s" % (compname, phasename)
                    ] = node_states.C[:, phase, comp]
                    dumped_states[
                        "fracture_comp_%s_in_phase_%s" % (compname, phasename)
                    ] = fracture_states.C[:, phase, comp]
        if dump_fluxes:
            cell_fluxes, fracture_fluxes = self.simulation.mass_fluxes()
            dumped_states["cell_total_mass_flux"] = cell_fluxes.sum(axis=1)
            dumped_states["fracture_total_mass_flux"] = fracture_fluxes.sum(axis=1)
            if True or nbcomponents > 1:
                for comp in range(nbcomponents):
                    compname = comp_names[comp]
                    dumped_states["cell_mass_flux_comp%s" % (compname)] = cell_fluxes[
                        :, comp, :
                    ]
                    dumped_states[
                        "fracture_mass_flux_comp%s" % (compname)
                    ] = fracture_fluxes[:, comp, :]
        np.savez(self.states_filename(mpi.proc_rank, tag), **dumped_states)
        dw.dump_all_wells(self.simulation, self.to_wells_directory(), tag)


def dump_mesh(simulation):
    dumper = Dumper(simulation)
    dumper.dump_own_element_numbers()
    dumper.dump_mesh()


def dump_states(simulation, tag=""):
    Dumper(simulation).dump_states(tag)
