# -*- coding: utf-8 -*-

# def mydecode(s, dtype='i'):
#     bs = base64.b64decode(s)
#     print np.fromstring(bs, dtype=dtype)

import sys, os
from functools import partial
import base64
import xml.dom
import numpy as np

vtk_celltype = ({
    'point': 1, 'vertex': 1, 'pt': 1,
    'line': 3,
    'triangle': 5, 'tri': 5,
    'pixel': 8,
    'quad': 9,
    'tetrahedron': 10, 'tet': 10,
    'voxel': 11,
    'hexahedron': 12, 'hex': 12,
    'wedge': 13,
})


def create_childnode(parent, name, attributes=None):
    if attributes is None:
        attributes = {}
    elt = parent.ownerDocument.createElement(name)
    for key, val in attributes.items():
        elt.setAttribute(key, val)
    parent.appendChild(elt)
    return elt


def vtk_doc(filetype):
    impl = xml.dom.getDOMImplementation()
    doc = impl.createDocument(None, 'VTKFile', None)
    doc.documentElement.setAttribute('type', filetype)
    doc.documentElement.setAttribute('header_type', 'UInt64')
    return doc

def add_xdataarray_node(parent, nodename, dataname, dtype, ofmt, nbcomp):
    attributes = {}
    if dtype == np.int32:
        dtype = 'Int32'
    elif dtype == np.int8:
        dtype = 'Int8'
    elif dtype == np.uint8:
        dtype = 'Int32'
    elif dtype == np.int64:
        dtype = 'Int64'
    elif dtype == np.uint32:
        dtype = 'Int32'
    elif dtype == np.uint64:
        dtype = 'Int64'
    elif dtype == np.float32:
        dtype = 'Float32'
    elif dtype == np.float64:
        dtype = 'Float64'
    elif dtype == np.double:
        dtype = 'Float64'
    else:
        raise Exception('Unknown data type: ' + str(dtype))
    attributes['Name'] = dataname
    attributes['type'] = dtype
    attributes['format'] = ofmt
    if nbcomp is not None:
        attributes['NumberOfComponents'] = '%d' % nbcomp
    return create_childnode(parent, nodename, attributes)

def add_dataarray_node(parent, name, dtype, ofmt='ascii', nbcomp=None):
    return add_xdataarray_node(parent, 'DataArray', name, dtype, ofmt, nbcomp)

def add_pdataarray_node(parent, name, dtype, ofmt='ascii', nbcomp=None):
    return add_xdataarray_node(parent, 'PDataArray', name, dtype, ofmt, nbcomp)

def add_dataarray(
    parent, array, name,
    ofmt='ascii', nbcomp=None,
    nbitemsbyrow=10,
):
    elt = add_dataarray_node(parent, name, array.dtype, ofmt, nbcomp)
    doc = elt.ownerDocument
    assert len(array.shape) == 1
    if ofmt == 'ascii':
        if array.dtype.kind in ['i', 'u']:
            datafmt = 'd'
        else:
            datafmt = '.10f'
        fmt = ' '.join(['{:%s}' % datafmt] * nbitemsbyrow)
        i = -1
        for i in range(0, array.shape[0] // int(nbitemsbyrow)):
            elt.appendChild(
                doc.createTextNode(
                    fmt.format(*array[i * nbitemsbyrow:(i + 1) * nbitemsbyrow])
                )
            )
        left = array[(i + 1) * nbitemsbyrow:]
        fmt = ' '.join(['{:%s}' % datafmt] * len(left))
        elt.appendChild(doc.createTextNode(fmt.format(*left)))
    elif ofmt == 'binary':
        elt.appendChild(doc.createTextNode(''))  # this is just to indent node
        nbytes = 8 + array.nbytes
        tmp = np.empty(nbytes, dtype=np.byte)
        tmp = np.require(tmp, requirements=['CONTIGUOUS', 'OWNDATA'])
        if sys.version_info.major<3:
            buffersize = np.empty(1, dtype=np.ulonglong)
            buffersize[0] = array.nbytes
            tmp[:8] = np.getbuffer(buffersize)
            tmp[8:] = np.getbuffer(array)
        else:
            buffersize = memoryview(tmp[:8]).cast('Q')
            buffersize[0] = array.nbytes
            datapart = memoryview(tmp[8:])
            datapart[:] = memoryview(array).cast('b')
        s = base64.b64encode(tmp).decode('ascii')
        elt.appendChild(doc.createTextNode(s))
    else:
        raise Exception('Unknown format !')
    return elt


def add_piece_data(piece, location, data=None, ofmt='ascii'):
    datanode = create_childnode(piece, location, {'Scalars': 'scalars'})
    if data is None:
        data = {}
    for name in sorted(data):
        a = data[name]
        assert len(a.shape)==1 or len(a.shape)==2
        if len(a.shape)==2:
            add_dataarray(datanode, _ravel_information_block(a), ofmt=ofmt, name=name, nbcomp=a.shape[1])
        else:
            add_dataarray(datanode, a, ofmt=ofmt, name=name)


def _ravel_information_block(block):
    return block.ravel(order='C')

# Paraview data array are in Fortran (column major) mode


def _ravel_structured_data_block(block):
    return block.ravel(order='F')


def vti_doc(
    shape, delta=None, origin=None,
    extent = None,
    pointdata=None, celldata=None,
    ofmt='binary',
):
    if extent is not None:
        assert origin is None and delta is None
        origin = tuple(x[0] for x in extent)
        # conversion to float in case the script is run with python2
        delta = tuple((x[1]-x[0])/float(nx) for nx, x in zip(shape, extent))
    if delta is None:
        assert extent is None
        delta = (1.,) * len(shape)
    assert len(shape) == len(delta)
    if origin is None:
        assert extent is None
        origin = (0.,) * len(shape)
    assert len(shape) == len(origin)
    if len(shape) < 3:
        nbfills = 3 - len(shape)
        shape = shape + (1,) * nbfills
        delta = delta + (1.,) * nbfills
        origin = origin + (0.,) * nbfills
    if pointdata is None:
        pointdata = {}
    else:
        pointdata = {key: _ravel_structured_data_block(block)
                     for key, block in pointdata.items()}
    if celldata is None:
        celldata = {}
    else:
        celldata = {key: _ravel_structured_data_block(block)
                    for key, block in celldata.items()}
    doc = vtk_doc('ImageData')
    s = ' '.join(['%f' for _ in origin])
    extent = ' '.join(['0 %d' for _ in shape]) % tuple(shape)
    grid = create_childnode(
        doc.documentElement, 'ImageData', {
            'Origin': s % tuple(origin),
            'Spacing': s % tuple(delta),
            'WholeExtent': extent
        })
    piece = create_childnode(grid, 'Piece', {'Extent': extent})
    add_piece_data(piece, 'PointData', pointdata, ofmt=ofmt)
    add_piece_data(piece, 'CellData', celldata, ofmt=ofmt)
    return doc

def vtu_vertices(vertices):
    vertices = np.array(vertices, copy=False)
    assert len(vertices.shape) == 2
    dim = vertices.shape[1]
    assert dim <= 3
    if dim < 3:
        tmp = np.zeros((vertices.shape[0], 3), dtype=vertices.dtype)
        tmp[:, :dim] = vertices
        vertices = tmp
    return vertices

def vtu_doc_from_COC(
    vertices, offsets, connectivity, celltypes,
    pointdata=None, celldata=None,
    ofmt='binary'
):
    offsets = offsets.astype(np.int64)
    celltypes = celltypes.astype(np.int8)
    doc = vtk_doc('UnstructuredGrid')
    grid = create_childnode(doc.documentElement, 'UnstructuredGrid')
    piece = create_childnode(
        grid, 'Piece', {
            'NumberOfPoints': '%d' % vertices.shape[0],
            'NumberOfCells': '%d' % celltypes.shape[0],
        },
    )
    points = create_childnode(piece, 'Points')
    add_dataarray(points, _ravel_information_block(vertices),
                  'Points', nbcomp=3, ofmt=ofmt)
    cells = create_childnode(piece, 'Cells')
    add_dataarray(cells,
                  _ravel_information_block(connectivity),
                  'connectivity', ofmt=ofmt)
    add_dataarray(cells, offsets, 'offsets', ofmt=ofmt)
    add_dataarray(cells, celltypes, 'types', ofmt=ofmt)
    add_piece_data(piece, 'PointData', pointdata, ofmt=ofmt)
    add_piece_data(piece, 'CellData', celldata, ofmt=ofmt)
    return doc

def vtu_doc(
    vertices, connectivity, celltypes=None,
    pointdata=None, celldata=None,
    ofmt='binary'
):
    vertices = vtu_vertices(vertices)
    try:
        connectivity = np.array(connectivity, dtype=np.int64, copy=False)
    except ValueError:
        # connectivities may have different length
        offsets = np.array(len(cell) for cell in connectivity)
        celltypes = np.array(celltypes, dtype='u1', copy=False)
    else:
        assert len(connectivity.shape) == 2
        nbcells, cellssize = connectivity.shape
        offsets = np.tile(cellssize, nbcells)
        if celltypes is None:
            celltypes = {
                1: 'pt', 2: 'line', 3: 'tri', 4: 'tet', 6: 'wedge', 8: 'hex',
            }[cellssize]
        if isinstance(celltypes, str):
            celltypes = np.tile(vtk_celltype[celltypes], nbcells)
    finally:
        offsets = offsets.astype(np.int64)
        offsets = np.cumsum(offsets)
        nbcells = offsets.shape[0]
    assert offsets.shape == celltypes.shape
    return vtu_doc_from_COC(
        vertices, offsets, connectivity, celltypes,
        pointdata, celldata,
        ofmt
    )

def pvtu_doc(
    vertices_type, pieces,
    pointdata_types=None, celldata_types=None,
    ofmt='binary'
):
    doc = vtk_doc('PUnstructuredGrid')
    root_element = doc.documentElement
    pugrid = create_childnode(root_element, 'PUnstructuredGrid', {'GhostLevel': '0'})
    ppoints = create_childnode(pugrid, 'PPoints')
    add_pdataarray_node(ppoints, 'Points', vertices_type, ofmt, nbcomp=3)
    def add_data_node(node_name, data_types):
        if data_types is not None:
            data_node = create_childnode(pugrid, node_name)
            for name in sorted(data_types):
                data_type = data_types[name]
                try:
                    data_type, nbcomp = data_type
                    add_pdataarray_node(data_node, name, data_type, ofmt, nbcomp=nbcomp)
                except TypeError:
                    add_pdataarray_node(data_node, name, data_type, ofmt)
    add_data_node('PPointData', pointdata_types)
    add_data_node('PCellData', celldata_types)
    for piece in pieces:
         create_childnode(pugrid, 'Piece', {'Source': piece})
    return doc

def points_as_vtu_doc(vertices, pointdata=None, ofmt='binary'):
    connectivity = np.arange(len(vertices))
    connectivity.shape = (-1, 1)
    return vtu_doc(
        vertices, connectivity, celltypes='point',
        pointdata=pointdata,
        ofmt=ofmt,
    )


def block_as_vti_doc(
    block, location, name,
    delta=None, origin=None,
    ofmt='binary',
):
    block = np.array(block, copy=False)
    pointdata = celldata = {}
    if location == 'point':
        pointdata = {name: _ravel_structured_data_block(block)}
        shape = tuple(n - 1 for n in block.shape)
    elif location == 'cell':
        celldata = {name: _ravel_structured_data_block(block)}
        shape = block.shape
    else:
        raise Exception('unknown location' + str(location))
    return vti_doc(
        shape, delta, origin, pointdata, celldata,
        ofmt=ofmt,
    )


def blocks_as_vti_doc(
    pointdata=None, celldata=None,
    delta=None, origin=None,
    ofmt='binary',
):
    if pointdata is None:
        pointdata = {}
    if celldata is None:
        celldata = {}
    pshapes = [v.shape for v in pointdata.values()]
    pshape = None if not pshapes else pshapes[0]
    if not(pshape is None or all([shape == pshape for shape in pshapes])):
        raise Exception('incompatible point data shapes')
    cshapes = [v.shape for v in celldata.values()]
    cshape = None if not cshapes else cshapes[0]
    if not(cshape is None or all([shape == cshape for shape in cshapes])):
        raise Exception('incompatible cell data shapes')
    if pshape is None and cshape is None:
        raise Exception('no data provided')
    if not(
        pshape is None or cshape is None or
        cshape == tuple(n - 1 for n in pshape)
    ):
        raise Exception('incompatible pointdata and celldata')
    if cshape is None:
        shape = tuple(n - 1 for n in pshape)
    else:
        shape = cshape
    return vti_doc(
        shape, pointdata=pointdata, celldata=celldata,
        delta=delta, origin=origin, ofmt=ofmt,
    )


def pvd_doc(snapshots):
    """ snapshots is assumed to be a collection of tuples with the format
    (time, filepath)"""
    doc = vtk_doc('Collection')
    collection = create_childnode(doc.documentElement, 'Collection')
    for t, filepath in sorted(snapshots):
        create_childnode(
            collection, 'DataSet',
            {'timestep': '%g' % t, 'file': str(filepath)},
        )
    return doc


def write_xml(doc, out, indent=' ' * 2, newl='\n', extension=''):
    def output(filelike):
        doc.writexml(filelike, addindent=indent, newl=newl)
    if isinstance(out, str):
        filename = out
        if not filename.endswith(extension):
            filename = out + extension
        with open(filename, 'w') as f:
            output(f)
        return filename
    else:
        output(out)
        return out


write_vti = partial(write_xml, extension='.vti')
write_vtu = partial(write_xml, extension='.vtu')
write_pvtu = partial(write_xml, extension='.pvtu')
write_pvd = partial(write_xml, extension='.pvd')


def _write_data_snapshots(
    doc_from_data,
    times, datas, name, filepath='.',
    propname=None, proppath='.', extension='',
    delta=None, origin=None, location='cell',
    ofmt='binary',
    indent=' ' * 2, newl='\n',
):
    if propname is None:
        propname = name
    filedir = os.path.join(filepath, proppath)
    if os.path.exists(filedir):
        assert os.path.isdir(filedir)
    else:
        os.makedirs(filedir)
    datapath = os.path.join(filedir, propname)
    datafiles = []
    # size of counter field
    fmt = '%%0%dd' % (int(np.log10(len(datas)) + 1))
    for i, data in enumerate(datas):
        datafiles.append(
            write_xml(
                doc_from_data(data, propname),
                datapath + (fmt % i),
                indent=indent, newl=newl,
                extension=extension,
            )
        )
    filepath = os.path.join(filepath, name)
    datafiles = [os.path.relpath(path, os.path.dirname(filepath))
                 for path in datafiles]
    return write_pvd(
        pvd_doc(list(zip(times, datafiles))),
        filepath, indent=indent, newl=newl,
    )


def write_block_snapshots(
    times, blocks, name, filepath='.',
    propname=None, proppath='.',
    delta=None, origin=None, location='cell',
    ofmt='binary',
    indent=' ' * 2, newl='\n',
):
    def doc_from_data(data, name):
        return block_as_vti_doc(
            data, location, name,
            delta=delta, origin=origin,
            ofmt=ofmt,
        )
    return _write_data_snapshots(
        doc_from_data, times, blocks, name,
        filepath=filepath, propname=propname, proppath=proppath,
        extension='.vti', indent=' ' * 2, newl='\n',
    )


def write_unstructured_snapshots(
    times, name, vertices, connectivity,
    datas, location,
    filepath='.', propname=None, proppath='.',
    celltypes=None, ofmt='binary',
    indent=' ' * 2, newl='\n',
):
    def doc_from_data(data, name):
        pointdata = None
        celldata = None
        datadict = {name: data}
        if location == 'point':
            pointdata = datadict
        if location == 'cell':
            celldata = datadict
        return vtu_doc(
            vertices, connectivity, celltypes=celltypes,
            pointdata=pointdata, celldata=celldata,
            ofmt=ofmt,
        )
    return _write_data_snapshots(
        doc_from_data, times, datas, name,
        filepath=filepath, propname=propname, proppath=proppath,
        extension='.vtu', indent=' ' * 2, newl='\n',
    )

def points_as_vtu(vertices):
    tmp = np.reshape(vertices, (-1, 3))
    return vtu_doc(
            tmp,
            np.rehsape(np.arange(tmp.shape[0]), (-1, 1)),
            )

if __name__ == '__main__':
    import itertools
    delta = 1., 1.5
    shape = 2, 3
    B = np.arange(np.prod(shape))
    B.shape = shape
    write_vti(block_as_vti_doc(B, location='cell', name='B'), 'image2D')
    delta = 1., 1.5, 2.
    shape = 2, 3, 4
    B1 = np.arange(np.prod(shape))
    B1.shape = shape
    write_vti(vti_doc(shape, delta, celldata={'B1': B1}), 'image1')
    write_vti(block_as_vti_doc(B1, location='point', name='B1'), 'image2')
    for shape in itertools.permutations((3, 5, 7)):
        size = np.prod(shape)
        block = np.zeros(shape, dtype='i')
        dim = len(shape)
        for axis in range(dim):
            axisname = ['Ox', 'Oy', 'Oz'][axis]
            filename = 'test-block_%dx%dx%d-along_' % (shape) + axisname
            pos = [slice(None), ] * dim
            for k in range(shape[axis]):
                pos[axis] = k
                block[tuple(pos)] = k
            doc = block_as_vti_doc(block, location='cell', name='foo')
            write_vti(doc, filename)
    B2 = 2 * B1
    doc = blocks_as_vti_doc(celldata={'B1': B1, 'B2': B2})
    write_vti(doc, 'image3')
    doc = blocks_as_vti_doc(
        pointdata={'B1': B1, 'B2': B2},
        celldata={'B1bis': B1[1:, 1:, 1:], 'B2bis': B2[1:, 1:, 1:]}
    )
    write_vti(doc, 'image4')
    # %% dummy 2D unstructured grid
    vertices = np.array([
        [[i, j] for i in range(3)] for j in range(3)
    ])
    vertices.shape = (-1, 2)
    connectivity = np.array([[0, 1, 3], [1, 3, 4], [1, 2, 4], [3, 4, 6]])
    write_vtu(vtu_doc(vertices, connectivity), 'grid')
    # %% block snapthots
    times = np.arange(10)
    blocks = [B1 + t for t in times]
    write_block_snapshots(times, blocks, 'foo', proppath='values')
    # %% unstructured snapshots
    C = np.arange(len(connectivity))
    datas = [C + t for t in times]
    write_unstructured_snapshots(
        times, 'ufoo', vertices, connectivity,
        datas, 'cell', proppath='values',
    )
