# -*- coding: utf-8 -*-
"""
Created on Tue Aug 04 14:33:26 2015

@author: lopez
"""

from itertools import combinations
import numpy as np

# this should be always the same as we use byte comparisons
_index_type = np.int

def _check_keyword(line, keyword, check=None, get_value_type=None):
    l = line.strip().split()
    if l[0].strip()!=keyword:
        raise Exception('Expecting keyword "%s"' % keyword)   
    if check is not None:
        if type(check)(l[1].strip())!=check:
            raise Exception('Expecting value %s for keyword "%s"'
                        % (str(check), keyword))   
    if get_value_type is not None:
        assert check is None
        try:
            value = get_value_type(l[1].strip())
        except ValueError:
            raise Exception('Could not read keyword "%s"' % (keyword))
        return value

def _read_indexed_table(f, keyword, dtype):
    _check_keyword(f.readline(), keyword)
    n = int(f.readline())
    data = np.array([f.readline().strip().split() for _ in range(n)],
                     dtype=dtype)
    index = np.array(data[:, -1], dtype='i', copy=True)
    data = np.array(data[:, :-1], copy=True)
    return data, index

def process(filename):    
    """Extract mesh information from a medit file.
    Vertices indexes are reset with 0 base instead of 1."""
    with open(filename) as f:
        version = _check_keyword(f.readline(), 'MeshVersionFormatted', get_value_type=int)
        if version==1:
            _check_keyword(f.readline(), 'Dimension', check=3)
        else:
            assert version==2
            _check_keyword(f.readline(), 'Dimension')
            dim = int(f.readline())
            assert dim==3
        vertices, vertices_index = _read_indexed_table(f, 'Vertices', np.double)    
        triangles, triangles_index = _read_indexed_table(f, 'Triangles', _index_type)    
        tets, tets_index = _read_indexed_table(f, 'Tetrahedra', _index_type)
        _check_keyword(f.readline(), 'End')
    # use index 0 for the first vertex instead of 1
    tets-= 1
    triangles-= 1
    assert tets.min()==0
    assert tets.max() + 1==vertices.shape[0]
    assert triangles.min()>=tets.min()
    assert triangles.max()<=tets.max()
    return (
        (vertices, vertices_index),
        (triangles, triangles_index),
        (tets, tets_index),
    )

def row_as_void_type(a):
    """abstract type used to compare rows as sequence of bytes"""
    assert len(a.shape)==2
    return np.dtype((np.void, a.dtype.itemsize * a.shape[1]))

def unique_rows(a, **kargs):
    """ return the index of a unique row of an array
        order of the array is supposed to be C order
        if the array is not contiguous a temporary copy is created
    """
    a = np.require(a, requirements='C')
    rowtype = row_as_void_type(a)
    #    res = np.unique(a.view(rowtype), return_index=True)
    kargscopy = dict(kargs)
    kargscopy['return_index'] = True    
    res = np.unique(a.view(rowtype), **kargscopy)
    idx = res[1]
    if 'return_index' in kargs and kargs['return_index']:
        return (a[idx],) + res[1:]
    return (a[idx],) + res[2:]

def rows_in1d(a, b):
    """ return the index of a unique row of an array
        order of the array is supposed to be C order
        if the array is not contiguous a temporary copy is created
    """
    assert len(a.shape)==len(b.shape)==2
    assert a.dtype==b.dtype
    a = np.require(a, 'C')
    b = np.require(b, 'C')
    rowtype = row_as_void_type(a)
    return np.in1d(
        a.view(rowtype), b.view(rowtype),
        assume_unique=assume_unique, invert=invert,
    )    

class MeditInfo(object):
    
    def __init__(self, filename):
        vertices_info, triangles_info, tets_info = process(filename)
        self.vertices, self.vertices_index = vertices_info        
        self.triangles, self.triangles_index = triangles_info        
        self.tets, self.tets_index = tets_info
        self.cells, self.cells_index = self.tets, self.tets_index

    def compute_face_info(self):
        cellfaces = np.vstack([
            np.array([face for face in combinations(tet, 3)], dtype=_index_type)
            for tet in self.tets            
        ])
        cellfaces.sort(axis=-1) # axis=-1 is the default
        faces, inverse = unique_rows(cellfaces, return_inverse=True)
        cellfaces = np.reshape(inverse, (-1, 4)) # 4 faces by tet
        facecells = [list() for _ in range(faces.shape[0])]
        for cj, cifaces in enumerate(cellfaces):
            for fi in cifaces:
                facecells[fi].append(cj)
        for twins in facecells:
            if len(twins)<2:
                assert len(twins)==1
                twins.append(-1)
        facecells = np.array(facecells, dtype=_index_type)
        self.faces = faces
        self.cellfaces = cellfaces
        self.facecells = facecells
        
    def faces_indexes(self, testfaces):
        """WARNING self.faces is assume to be unique"""
        faces = self.faces
        testfaces = np.array(testfaces, dtype=_index_type, copy=True)
        testfaces.shape = (-1, 3)        
        testfaces.sort(axis=-1)
        allfaces = np.vstack([faces, testfaces])
        _, idx = unique_rows(allfaces, return_inverse=True)
        nbfaces = faces.shape[0] # this is the maximum facet index
        testidx = idx[nbfaces:]
        testidx[testidx>=nbfaces] = -1
        return testidx

    def compute_edge_info(self):
        try:
            faces = self.faces
        except AttributeError:
            self.compute_face_info()
            faces = self.faces
        faceedges = np.vstack([
            np.array([edge for edge in combinations(face, 2)], dtype=_index_type)
            for face in faces            
        ])
        faceedges.sort(axis=-1) # axis=-1 is the default
        edges, inverse = unique_rows(faceedges, return_inverse=True)
        faceedges = np.reshape(inverse, (-1, 3)) # 3 edges by face
        self.edges = np.array([list(edge) for edge in edges])
        self.faceedges = faceedges
    
    def compute_info(self):
        self.compute_face_info()
        self.compute_edge_info()

#    def output_vtu(self, basename, output_triangles=True):
#        vtkw.write_vtu(
#            vtkw.vtu_doc(self.vertices, self.cells,
#                         celldata={'material': self.cells_index}),
#            basename,
#        )
#        if output_triangles:
#            absolute_triangle_index = self.faces_indexes(self.triangles)
#            bad_triangles = np.array(absolute_triangle_index==-1, dtype=np.int)     
#            vtkw.write_vtu(
#                vtkw.vtu_doc(self.vertices, self.triangles,
#                             celldata={'surface_patch': self.triangles_index,
#                                       'bad': bad_triangles}),
#                basename + '_triangles',
#            )

def reader(filename):
    return MeditInfo(filename)
     
