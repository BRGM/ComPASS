# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 08:51:49 2018

@author: lopez
"""

import sys
import os
import re
import numpy as np
from scipy.sparse import csr_matrix
import scipy.sparse.linalg as spl

# to control the use of atol for gmres
import scipy.version as scpver
scipy_version = int(scpver.full_version.split('.')[1])

if len(sys.argv)!=2:
    print('to few arguments')
    sys.exit(-1)

basename = sys.argv[1]

Afile = '%s_A.dat' % basename
if not os.path.exists(Afile):
    print(Afile, 'not found')
    sys.exit(-1)

bfile = '%s_b.dat' % basename
if not os.path.exists(bfile):
    print(bfile, 'not found')
    sys.exit(-1)

filter_out_zeros = True

matcol = []
matdat = []
with open(Afile) as f:
    line = f.readline().strip()
    match = re.match('Mat Object: (\d+) MPI processes', line)
    # 1 proc only
    assert match and int(match.groups()[0])==1
    f.readline()
    column_extraction = re.compile('row (\d+): (.*)')
    data_extraction = re.compile('\((\d+), (.+)\)')
    line = f.readline().strip()
    while line:
        match = column_extraction.match(line)
        groups = match.groups()
        row = int(groups[0])
        data = [s.strip() for s in groups[1].strip().split('  ')]
        cols = []
        nz = []
        for s in data:
            groups = data_extraction.match(s).groups()
            cols.append(int(groups[0]))
            try:
                nz.append(float(groups[1]))
            except ValueError:
                print('\nconversion error on row %d:\n%s\n' % (row, data))
                raise
        assert len(matcol) == len(matdat) == row
        cols = np.array(cols)
        nz = np.array(nz)
        assert not np.any(np.isnan(nz))
        if filter_out_zeros:
            keep = nz!=0
        else:
            keep = np.ones(nz.shape[0], dtype='b')
        matcol.append(cols[keep])
        matdat.append(np.array(nz[keep]))
        line = f.readline().strip()

lengths = np.array([a.shape[0] for a in matcol])
assert np.all(np.array([a.shape[0] for a in matdat])==lengths)

ptr = np.hstack([[0], np.cumsum(lengths)])

assert len(matdat)==row+1
A = csr_matrix((
        np.hstack(matdat), np.hstack(matcol),
        np.hstack([[0], np.cumsum(lengths)]),
        ), shape=(row+1, row+1))
print('matrix size', A.shape)

bvalues = []
with open(bfile) as f:
    for _ in range(3): # header
        f.readline()
    line = f.readline().strip()
    while line:
        bvalues.append(float(line))
        line = f.readline().strip()
b = np.array(bvalues)
assert b.shape[0]==row+1

xdirect = spl.spsolve(A, b)

norm = np.linalg.norm

print('Direct solution ', norm(A.dot(xdirect)-b))
assert np.allclose(A.dot(xdirect), b)

residuals_gmres = []
tol, atol = 1e-8, 1e-10
maxiter = 200
if scipy_version<20:
    xgmres, return_code = spl.gmres(
        A, b, tol=tol, restart=20, maxiter=maxiter,
        callback=lambda rk: residuals_gmres.append(rk)
    )
else:
    xgmres, return_code = spl.gmres(
        A, b, tol=tol, atol=atol, restart=20, maxiter=maxiter,
        callback=lambda rk: residuals_gmres.append(rk)
    )
if return_code==0: #convergence
    assert norm(A.dot(xgmres) - b)<max(tol*norm(b), atol)
    print('Convergence gmres in', len(residuals_gmres), 'iterations')
    print('GMRES solution ', norm(A.dot(xgmres)-b))
else:
    print('GMRES did not converged!')
    if return_code>0:
        print('Iteration exhaustion:', return_code)
#print('Solutions differ by', norm(xdirect-xgmres))

#residuals_bicgstab = []
#xbicgstab, return_code = spl.gmres(
#        A, b, tol=tol,
#        callback=lambda rk: residuals_bicgstab.append(rk)
#)
#assert return_code==0 #convergence
#print('Convergence bicgstab in', len(residuals_bicgstab), 'iterations')
#print('BICGSTAB solution ', norm(A.dot(xbicgstab)-b))
#print('Solutions differ by', norm(xdirect-xbicgstab))

#import matplotlib.pyplot as plt
#plt.plot(np.arange(len(residuals_gmres)), residuals_gmres)
#plt.plot(np.arange(len(residuals_bicgstab)), residuals_bicgstab)


norm_A = spl.norm(A)
norm_invA = spl.norm(spl.inv(A.tocsc()))
cond = norm_A*norm_invA
print('condition number %.2e' % cond)

def explore(s, v):
    av = np.abs(v)
    print(s, np.min(av), np.mean(av), np.max(av))
explore('b p part', b[::2])
explore('b T part', b[1::2])

