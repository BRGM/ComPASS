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
import numpy as np
import thirdparties
import vtkwriters as vtkw

def collect_snapshots(snapshots_file):
    assert os.path.exists(snapshots_file)
    with open(snapshots_file) as f:
        snapshots = {}
        for line in f:

            l = line.strip().split()
            if len(l)>=2:
                snapshots[l[0]] = float(l[1])
    return snapshots

def collect_procs(owns_file):
    assert os.path.exists(owns_file)
    with open(owns_file) as f:
        line = f.readline()
        while line and not line.strip().lower().startswith('nb_procs ='):
            line = f.readline()
        line = line.strip().lower()
        result = re.match('nb_procs =\s*(\d*)', line)
        assert result
        return int(result.group(1))

def create_directories(path):
    if not os.path.exists(path):
        os.makedirs(path)
    assert os.path.isdir(path)

def extract_data(datafile):
    assert os.path.exists(datafile)
    data = np.load(datafile)
    def extract_at(location):
        return { 
                  name.replace(location + '_', ''): data[name]
                  for name in data.files if name.startswith(location)
               }
    return {
              location: extract_at(location)
              for location in ('node', 'cell', 'fracture')
           }

def collect(directory):
    assert os.path.isdir(directory)
    print('processing results in', directory)
    todir = lambda s='' : os.path.join(directory, s)
    topardir = lambda s='' : os.path.join(directory, 'paraview', s)
    create_directories(topardir())
    toparvtudir = lambda s='' : os.path.join(topardir('vtus'), s)
    create_directories(toparvtudir())
    snapshots = collect_snapshots(todir('snapshots'))
    nbprocs = collect_procs(todir('own_elements'))
    pvdcollection = []
    fracpvdcollection = []
    for tag in snapshots:
        #print('processing tag', tag)
        vertices_type = None
        data_types = None
        pieces = []
        fracpieces = []
        for p in range(nbprocs):
            proc = 'proc_%04d' % p
            meshfile = todir(os.path.join('mesh', 'mesh_' + proc + '.npz'))
            assert os.path.exists(meshfile)
            mesh = np.load(meshfile)
            datafile_basename = 'state_' + tag +'_' + proc + '.npz'
            datafile = todir(os.path.join('states', 'output-' + tag, datafile_basename))
            data_arrays = extract_data(datafile)
            datatype_checks = {
                                location: {
                                        name: property.dtype
                                        for name, property in properties.items()
                                    } for location, properties in data_arrays.items()
                              }
            if vertices_type is None:
                vertices_type = mesh['vertices'].dtype
            assert vertices_type == mesh['vertices'].dtype
            if data_types is None:
                data_types = datatype_checks
            assert data_types == datatype_checks
            piecefile = datafile_basename.replace('.npz', '.vtu')
            vtkw.write_vtu(
                vtkw.vtu_doc_from_COC(
                    mesh['vertices'], 
                    mesh['cellnodes_offsets'], 
                    mesh['cellnodes_values'], 
                    mesh['celltypes'], 
                    pointdata = data_arrays['node'],
                    celldata = data_arrays['cell'],
                ),
                toparvtudir(piecefile)
            )
            pieces.append(piecefile)
            if len(mesh['fracture_types'])>0:
                fracpiecefile = 'fracture_' + piecefile
                vtkw.write_vtu(
                    vtkw.vtu_doc_from_COC(
                        mesh['vertices'], 
                        mesh['fracturenodes_offsets'], 
                        mesh['fracturenodes_values'], 
                        mesh['fracture_types'], 
                        celldata = data_arrays['fracture'],
                    ),
                    toparvtudir(fracpiecefile)
                )
                fracpieces.append(fracpiecefile)
        pvtufile = toparvtudir('state_' + tag) + '.pvtu'
        pvdcollection.append((snapshots[tag] / (86400*365.25), # one year
                              os.path.relpath(pvtufile, start=topardir())))
        vtkw.write_pvtu(
            vtkw.pvtu_doc(
                vertices_type, pieces,
                pointdata_types = data_types['node'],
                celldata_types = data_types['cell'],
            ),
            pvtufile
        )
        if fracpieces:
            fracpvtufile = toparvtudir('fracture_' + 'state_' + tag) + '.pvtu'
            fracpvdcollection.append((snapshots[tag] / (86400*365.25), # one year
                                  os.path.relpath(fracpvtufile, start=topardir())))
            vtkw.write_pvtu(
                vtkw.pvtu_doc(
                    vertices_type, fracpieces,
                    celldata_types = data_types['fracture'],
                ),
                fracpvtufile
            )
    vtkw.write_pvd(
        vtkw.pvd_doc(pvdcollection),
        topardir('mesh_states')
    )
    if fracpvdcollection:
        vtkw.write_pvd(
            vtkw.pvd_doc(fracpvdcollection),
            topardir('fracture_states')
        )
    

if len(sys.argv)>1:
    for directory in sys.argv[1:]:
        collect(directory)
