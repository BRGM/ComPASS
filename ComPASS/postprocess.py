# -*- coding: utf-8 -*-
#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import sys
import glob, os, re
from optparse import OptionParser
import numpy as np
from .dumps import Dumper
from .utils import create_directories
from .utils.units import year
import MeshTools.vtkwriters as vtkw


parser = OptionParser(
usage = """usage: %prog [options] path

NB: you may also use this script as a module from the python package, i.e. typically:
python -m ComPASS.postprocesss [options] path

Parameters:
  path           where output from a ComPASS simulation are stored (typically output-scriptname)"""
)
parser.add_option("-p", "--procs",
                  action="store_true", dest="collect_procs_id", default=False,
                  help="ouput paraview/mesh.pvtu file with cell distribution (proc variable)")
parser.add_option("-s", "--states",
                  action="store_true", dest="collect_states", default=False,
                  help="ouput paraview files with transient states as pvd files to be found in the path/paraview directory")
parser.add_option("-C", "--Celsius",
                  action="store_true", dest="convert_temperature", default=False,
                  help="convert temperature from Kelvin to Celsius degrees")
options, args = parser.parse_args()

class MeshDistribution:
    def __init__(self, filename):
        def extract_integer_list(s):
            match = re.match('\[(.*)\]', s.strip())
            assert match is not None
            return  [int(si.strip()) for si in match.group(1).split()]
        vars = {
            'nb_procs': lambda s: int(s),
            'nb_own_cells': extract_integer_list,
            'nb_own_nodes': extract_integer_list,
            'nb_own_faces': extract_integer_list,
            'nb_own_fractures': extract_integer_list,
        }
        for name in vars:
            setattr(self, name, None)
        with open(filename) as f:
            for line in f:
                l = line.lower().strip()
                for name, extract in vars.items():
                    if getattr(self, name) is None:
                        match = re.match(name + '\s*=\s*(.*)', l)
                        if match is not None:
                            setattr(self, name, extract(match.group(1)))
                            break
    def __str__(self):
        return '\n'.join(
            [
                '%s = %d'%(s, getattr(self, s)) for s in (
                    nb_procs, nb_own_cells, nb_own_nodes,
                    nb_own_faces, nb_own_fractures
                )
            ]
        )

class PostProcessor:

    def __init__(self, directory):
        self.dumper = Dumper(directory)
        self.paraview_directory = os.path.join(
                                 self.dumper.to_output_directory(), 'paraview')
        create_directories(self.paraview_directory)
        self.vtu_directory = os.path.join(self.paraview_directory, 'vtu')
        create_directories(self.vtu_directory)
        self.distribution = MeshDistribution(self.dumper.to_output_directory('own_elements'))
        self.vertices_type = None
        self.data_types = None
        self.proc_id_type = np.dtype('i')

    def to_paraview_directory(self, filename=''):
        return os.path.join(self.paraview_directory, filename)
        
    def to_vtu_directory(self, filename=''):
        return os.path.join(self.vtu_directory, filename)
    
    def mesh_filename(self, proc):
        assert proc < self.distribution.nb_procs
        return self.dumper.mesh_filename(proc)
        
    def states_filename(self, proc, tag=''):
        assert proc < self.distribution.nb_procs
        return self.dumper.states_filename(proc, tag)
        
    def write_mesh_vtu(self, proc, basename,
                       own_only=True, dump_procid=False,
                       nodedata=None, celldata=None, fracdata=None):
        assert proc < self.distribution.nb_procs
        nb_own_cells = self.distribution.nb_own_cells[proc]
        nb_own_fractures = self.distribution.nb_own_fractures[proc]
        if dump_procid:
            own_only = True
            assert nodedata is None and celldata is None and fracdata is None
            celldata = {'proc': np.array(np.tile(proc, nb_own_cells),
                                         dtype=self.proc_id_type)}
            if nb_own_fractures>0:
                fracdata = {'proc': np.array(np.tile(proc, nb_own_fractures),
                                             dtype=self.proc_id_type)}
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
            self.vertices_type = mesh['vertices'].dtype
        assert self.vertices_type == mesh['vertices'].dtype
        proc_label = self.dumper.proc_label(proc)
        piecefile = self.to_vtu_directory('%s_%s.vtu' % (basename, proc_label))
        cell_nodes_offsets = mesh['cellnodes_offsets']
        cell_nodes = mesh['cellnodes_values']
        cell_types = mesh['celltypes']
        if own_only:
            cell_nodes_offsets = cell_nodes_offsets[:nb_own_cells]
            cell_nodes = cell_nodes[:cell_nodes_offsets[-1]]
            cell_types = cell_types[:nb_own_cells]
        vtkw.write_vtu(
            vtkw.vtu_doc_from_COC(
                mesh['vertices'], cell_nodes_offsets, cell_nodes, cell_types, 
                pointdata = nodedata, celldata=celldata
            ),
            piecefile
        )
        if fracdata:
            fracdata_size = int(np.unique([a.shape[0] for a in fracdata.values()])) # will triger a TypeError if array sizes are not the same
            if fracdata_size>0:
                fracpiecefile = self.to_vtu_directory('fracture_%s_%s.vtu' % (basename, proc_label))
                cell_nodes_offsets = mesh['fracturenodes_offsets']
                cell_nodes = mesh['fracturenodes_values']
                cell_types = mesh['fracture_types']
                if own_only:
                    cell_nodes_offsets = cell_nodes_offsets[:nb_own_fractures]
                    cell_nodes = cell_nodes[:cell_nodes_offsets[-1]]
                    cell_types = cell_types[:nb_own_fractures]
                vtkw.write_vtu(
                    vtkw.vtu_doc_from_COC(
                        mesh['vertices'], cell_nodes_offsets, cell_nodes, cell_types, 
                        celldata=fracdata
                    ),
                    fracpiecefile,
                )
                return piecefile, fracpiecefile
        return piecefile, None
    
    def collect_proc_ids(self):
        pieces, fracpieces = [], []
        for proc in range(self.distribution.nb_procs):
            piecefiles = self.write_mesh_vtu(proc, 'mesh', dump_procid=True)
            assert piecefiles[0] is not None
            pieces.append(os.path.relpath(piecefiles[0], self.to_paraview_directory()))
            if piecefiles[1] is not None:
                fracpieces.append(os.path.relpath(piecefiles[1], self.to_paraview_directory()))
        vtkw.write_pvtu(
            vtkw.pvtu_doc(
                self.vertices_type, pieces,
                celldata_types = {'proc': self.proc_id_type},
            ),
            self.to_paraview_directory('mesh.pvtu')
        )
        if len(fracpieces)>0:
            vtkw.write_pvtu(
                vtkw.pvtu_doc(
                    self.vertices_type, fracpieces,
                    celldata_types = {'proc': self.proc_id_type},
                ),
                self.to_paraview_directory('fractures_mesh.pvtu')
            )  
        
    def collect_snapshots(self):
        snapshots_file = self.dumper.to_output_directory('snapshots')
        assert os.path.exists(snapshots_file)
        with open(snapshots_file) as f:
            snapshots = {}
            for line in f:
    
                l = line.strip().split()
                if len(l)>=2:
                    snapshots[l[0]] = float(l[1])
        return snapshots

    def extract_data(self, proc, tag):
        datafile = self.states_filename(proc, tag)
        assert os.path.exists(datafile)
        data = np.load(datafile)
        def extract_at(location):
            return { 
                      name.replace(location + '_', ''): data[name]
                      for name in data.files if name.startswith(location)
                   }
        data_arrays = {
                  location: extract_at(location)
                  for location in ('node', 'cell', 'fracture')
               }
        datatype_checks = {
                            location: {
                                    name: property.dtype if len(property.shape)==1 else (property.dtype, property.shape[1])
                                    for name, property in properties.items()
                                } for location, properties in data_arrays.items()
                          }
        if self.data_types is None:
            self.data_types = datatype_checks
        assert self.data_types == datatype_checks
        return data_arrays
    
    def collect_states(self, options):
        snapshots = self.collect_snapshots()
        pvd = {}
        for tag in snapshots:
            all_pieces = []
            for proc in range(self.distribution.nb_procs):
                data_arrays = self.extract_data(proc, tag)
                if options.convert_temperature:
                    for location in data_arrays.keys():
                        for property in data_arrays[location].keys():
                            if property=='temperature':
                                data_arrays[location][property]-= 273.15
                pieces = self.write_mesh_vtu(proc, 'states_' + tag,
                                             nodedata=data_arrays['node'],
                                             celldata=data_arrays['cell'],
                                             fracdata=data_arrays['fracture'])
                all_pieces.append(pieces)
            pvtufile = self.to_vtu_directory('state_' + tag + '.pvtu')            
            vtkw.write_pvtu(
                vtkw.pvtu_doc(
                    self.vertices_type,
                    [os.path.relpath(pieces[0], self.to_vtu_directory()) for pieces in all_pieces],
                    pointdata_types = self.data_types['node'],
                    celldata_types = self.data_types['cell']
                ),
                pvtufile
            )
            pvtus = (pvtufile,)
            if any(pieces[1] is not None for pieces in all_pieces):
                fracpvtufile = self.to_vtu_directory('fracture_state_' + tag + '.pvtu')            
                vtkw.write_pvtu(
                    vtkw.pvtu_doc(
                        self.vertices_type,
                        [os.path.relpath(pieces[1], self.to_vtu_directory()) for pieces in all_pieces if pieces[1] is not None],
                        celldata_types = self.data_types['fracture']
                    ),
                    fracpvtufile
                )
                pvtus = (pvtufile, fracpvtufile)
            pvd[snapshots[tag] / year] = pvtus
        vtkw.write_pvd(
            vtkw.pvd_doc(
                [(t, os.path.relpath(pvtus[0], self.to_paraview_directory()))
                    for t, pvtus in pvd.items()]                    
            ),
            self.to_paraview_directory('states')
        )
        if any(len(pvtus)>1 for pvtus in pvd.values()):
            assert all(len(pvtus)>1 for pvtus in pvd.values())
            vtkw.write_pvd(
                vtkw.pvd_doc(
                    [(t, os.path.relpath(pvtus[1], self.to_paraview_directory()))
                        for t, pvtus in pvd.items()]                    
                ),
                self.to_paraview_directory('fracture_states')
            )

for directory in args:
    print('processing results in', directory)
    something_done = False
    if options.collect_procs_id or options.collect_states:
        pp = PostProcessor(directory)
        if options.collect_procs_id:
            pp.collect_proc_ids()
            something_done = True
        if options.collect_states:
            pp.collect_states(options)
            something_done = True
    if not something_done:
        print()
        print('************* WARNING **************')
        print()
        print('Nothing done: ckeck command options!')
        print()
        print('************************************')
        print()
