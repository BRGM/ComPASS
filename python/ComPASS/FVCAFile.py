import os
import numpy as np


class FVCAFile:

    def __init__(self, filename):
        assert os.path.exists(filename)
        with open(filename) as f:
            def skip(n):
                for _ in range(n):
                    f.readline()
            def extract_tag():
                f.readline().strip().lower()
            def read_and_check_tag(tag):
                return extract_tag()==tag
            def extract_connectivity():
                convert_strings = lambda ss: [int(s)-1 for s in ss] # -1: Fortran to C indexing
                l = f.readline().strip().split()
                n = int(l[0])
                result = convert_strings(l[1:])
                n-= len(result)
                while n>0:
                    l = f.readline().strip().split()
                    n-= len(l)
                    result.extend(convert_strings(l))
                assert n==0
                return result
            def tagged_int():
                l = f.readline().strip().split()
                tag = " ".join(l[:-1])
                n = int(l[-1])
                print('*', n, tag)
                return tag, n
            def find_tagged_int(tag):
                tag = tag.lower()
                l = []
                while len(l)!=2 and not l[0].strip().lower().startswith(tag):
                    l = f.readline().strip().split()
                tag, n = l[0].strip(), int(l[1])
                return tag, n
            # ---
            skip(9)
            nbnodes = int(f.readline().strip())
            skip(1)
            nbcells = int(f.readline().strip())
            skip(1)
            nbfaces = int(f.readline().strip())
            skip(1)
            nbedges = int(f.readline().strip())
            _, n = tagged_int()
            assert n==nbnodes
            self.vertices = np.array([
                    [float(s) for s in f.readline().split()]
                    for _ in range(nbnodes)
                ], dtype=np.double)
            _, n = tagged_int()
            assert n==nbcells
            self.cell_faces = [extract_connectivity() for _ in range(nbcells)]
            _, n = tagged_int()
            assert n==nbcells
            self.cell_nodes = [extract_connectivity() for _ in range(nbcells)]
            _, n = tagged_int()
            assert n==nbfaces
            self.face_edges = [extract_connectivity() for _ in range(nbfaces)]
            _, n = tagged_int()
            assert n==nbfaces
            self.face_nodes = [extract_connectivity() for _ in range(nbfaces)]
            _, n = tagged_int()
            assert n==nbfaces
            self.face_cells = [[int(s) for s in f.readline().strip().split()] for _ in range(nbfaces)]
            _, n = tagged_int()
            assert n==nbedges
            self.edge_nodes = np.array([
                    [int(s) for s in f.readline().strip().split()]
                    for _ in range(nbedges)
                ])
