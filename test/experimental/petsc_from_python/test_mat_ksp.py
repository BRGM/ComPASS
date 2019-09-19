import ComPASS.utils.mpl_backends as mplb

plt = mplb.import_pyplot()

try:
    execfile
except NameError:
    def execfile(file, globals=globals(), locals=locals()):
        fh = open(file, "r")
        try: exec(fh.read()+"\n", globals, locals)
        finally: fh.close()

import petsc4py, sys
petsc4py.init(sys.argv)

from petsc4py import PETSc

try: range = xrange
except: pass

from petsc4py import PETSc

# grid size and spacing
m, n  = 32, 32
hx = 1.0/(m-1)
hy = 1.0/(n-1)

# create sparse matrix
A = PETSc.Mat()
A.create(PETSc.COMM_WORLD)
A.setSizes([m*n, m*n])
A.setType('aij') # sparse
A.setPreallocationNNZ(5)

# precompute values for setting
# diagonal and non-diagonal entries
diagv = 2.0/hx**2 + 2.0/hy**2
offdx = -1.0/hx**2
offdy = -1.0/hy**2

# loop over owned block of rows on this
# processor and insert entry values
Istart, Iend = A.getOwnershipRange()
for I in range(Istart, Iend) :
    A[I,I] = diagv
    i = I//n    # map row number to
    j = I - i*n # grid coordinates
    if i> 0  : J = I-n; A[I,J] = offdx
    if i< m-1: J = I+n; A[I,J] = offdx
    if j> 0  : J = I-1; A[I,J] = offdy
    if j< n-1: J = I+1; A[I,J] = offdy

# communicate off-processor values
# and setup internal data structures
# for performing parallel operations
A.assemblyBegin()
A.assemblyEnd()


from ComPASS.bind_petsc import dump
dump(A)

# create linear solver
ksp = PETSc.KSP()
ksp.create(PETSc.COMM_WORLD)
# use conjugate gradients
ksp.setType('cg')
# and incomplete Cholesky
ksp.getPC().setType('icc')
# obtain sol & rhs vectors
x, b = A.createVecs()
x.set(0)
b.set(1)
# and next solve
ksp.setOperators(A)
ksp.setFromOptions()
ksp.solve(b, x)

OptDB = PETSc.Options()

if OptDB.getBool('plot_mpl', False):
    from numpy import mgrid
    X, Y =  mgrid[0:1:1j*m,0:1:1j*n]
    Z = x[...].reshape(m,n)
    plt.figure()
    plt.contourf(X,Y,Z)
    plt.plot(X.ravel(),Y.ravel(),'.k')
    plt.axis('equal')
    plt.colorbar()
    plt.savefig(__file__+'.png')
