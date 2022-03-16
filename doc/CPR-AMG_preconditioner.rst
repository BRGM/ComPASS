.. _cpramg:

The CPR-AMG preconditioner
=======================================

Introduction
-----------------

Linear systems in ComPASS describe complex physics on mesh degrees of
freedom of different natures (reservoir pressure or temperature,
wells head pressure etc..), governed
by partial differential equations with distinct analytical properties.
Usual preconditioning methods like ILU, Schwarz methods and multigrid
algorithms prove to be unsufficient on large systems because they
make no distinction between the unknowns,
thus losing knowledge of their physical and
numerical properties. The default preconditioner in ComPASS is CPR-AMG,
a method which uses the splitting of the blocks of unknowns in the linear
operator to improve the convergence of Krylov methods. This document
describes the CPR-AMG preconditioning procedure and its Python
implementation in the linear algebra package for
ComPASS.

The CPR-AMG preconditioning method
----------------------------------

Let :math:`A` be a linear operator of size :math:`n \times n` describing
the behaviour of two physical quantities, the pressure :math:`p`
and temperature :math:`T` (the matter here is to split the unknowns into a pressure block and
the rest of the unknowns, :math:`T` being a representant of the possible variables),
such that :math:`A` writes :

.. math::

   A=\left(\begin{array}{ll}
   A_{p p} & A_{pT} \\
   A_{T p} & A_{T T}
   \end{array}\right).

The aim is to solve the linear system :math:`A x=b`
using the CPR (Constrained Pressure Residual) preconditioner
([Roy_2019]_).
This is a left preconditioner, thus the preconditioned system
writes :math:`P^{-1} A x = P^{-1} b`. Matrix :math:`P^{-1}` aims at being an approximation
of :math:`A^{-1}` to reduce the system's condition number which
improves the convergence of Krylov methods. In practice a preconditioner can be entirely
defined by its application to a vector `r`, e.g. the value of :math:`u = P^{-1} r`.

CPR application is performed using two intermediate vectors :math:`v` and
:math:`w` of size :math:`n` as follows:

.. math::

      \begin{array}{l}
      v = 0, \\
      v_{|p} = \tilde{A}_{p p}^{-1} r_{|p}, \\
      w = r - A v, \\
      u = v + \tilde{A}^{-1} w,
      \end{array}

where :math:`v_{|p}` (respectively :math:`r_{|p}`) is the restriction vector
to the pressure unknowns of :math:`v` (resp. :math:`r`) ie
:math:`v_{|p} = R_p v` (resp. :math:`r_{|p} = R_p r`), with :math:`R_p`
the restriction matrix to the pressure unknowns.
(such that :math:`R_p A R_p^T = A_{p p}`). The approximation of :math:`A_{p p}^{-1}`,
denoted :math:`\tilde{A}_{p p}^{-1}`, is obtained by an algebraic multigrid
procedure and :math:`\tilde{A}^{-1}` approaches :math:`A^{-1}` using an ILU
pseudo-solve. This algorithm corresponds to multiplying :math:`r`
by the following matrix

.. math::

   \text{CPR-AMG}(A) = P^{-1} = M_2^{-1}\left(I-A M_1^{-1}\right) + M_1^{-1}

with the two sub-preconditioners
:math:`M_1^{-1}=\left(\begin{array}{cc} \tilde{A}_{p p}^{-1} & 0 \\ 0 & 0 \end{array}\right)` and
:math:`M_2^{-1} = \tilde{A}^{-1}`.

Implementation with petsc4py
----------------------------

Preconditioners of the form :math:`P^{-1} = M_2^{-1}\left(I-A M_1^{-1}\right)  + M_1^{-1}`
like CPR-AMG are an extension of the
Gauss-Seidel preconditioners where :math:`M_1^{-1}` and :math:`M_2^{-1}` can be any type of
preconditioning matrices. In PETSc they are implemented under the type ``COMPOSITE``
in ``MULTIPLICATIVE`` mode (additive application corresponds to the matrix
:math:`P^{-1} = M_2^{-1} + M_1^{-1}`).
The corresponding documentation can be found `in the PETSc manual <https://petsc.org/release/docs/manual/ksp/?highlight=preconditioner#combining-preconditioners>`_.
The PETSc implementation of the CPR-AMG preconditioning procedure
is done in file :download:`preconditioners.py <../ComPASS/linalg/preconditioners.py>`, follows
the corresponding extract.

.. literalinclude:: ../ComPASS/linalg/preconditioners.py
  :language: python
  :pyobject: CPRAMG
  :start-at: class CPRAMG(PETSc.PC):
  :end-at: self.setUp()
  :caption: Composite structure of CPR-AMG and its inner preconditioners type set with the Python interface to PETSc

One then has to set the internal sub-preconditioners
:math:`M_1^{-1}` and :math:`M_2^{-1}` types. Computing :math:`M_2^{-1}` is very straight forward because it is
a PETSc built-in block Jacobi algorithm which consists in applying
an incomplete LU solve on the local diagonal matrix block of each process. It is highly
scalable since entirely local, but the quality of the preconditioner deteriorates when the
number of processors increases.
The AMG procedure :math:`M_1^{-1}` aims at counteracting this effect by using knowledge
of the pressure equation system on the whole domain. Pseudo code for applying :math:`M_1^{-1}`
to a vector :math:`r` can be found below.

.. code::
  :caption: Application of the multigrid preconditioner on the pressure unknowns

  function apply_M1(A, r):
    set vector u to zero
    rp = Rp*r
    Ap <= Rp*A*RpT
    Mp <= amg_procedure(Ap)
    pressure restriction of u <= Mp*rp
    return u

The petsc4py API allows the user to set its own preconditioner using a Python class
defining the way the preconditioner must be set up and applied to a vector.
In this case the preconditioner must use the restriction of operator A to the
pressure equations. It is then necessary to identify a set of indexes,
it is done at ``setUp`` stage using a PETSc IndexSet object
(see `PETSc user's guide <https://petsc.org/release/docs/manual/vec/#index-sets>`_) for further details).

The ``apply`` method then retrieves the restriction of the pressure unknowns of the input
vector and copies it in the pressure part of the output vector.
Other indices of the output vector are null.

.. literalinclude:: ../ComPASS/linalg/preconditioners.py
  :language: python
  :pyobject: PressureAMG
  :start-at: def setUp(self, pc):
  :end-at: y.restoreSubVector(self.pressure_IS, return_subvec)



References
-----------------

.. [Roy_2019]  *A block preconditioner for non-isothermal flow in porous media*;
                 Thomas Roy, Tom B. Jönsthövel, Christopher Lemon and Andrew J. Wathen;
                 October 2019
                 http://dx.doi.org/10.1016/j.jcp.2019.06.038

.. [PETSc_manual] *PETSc Users Manual*; Satish Balay et al.;
                   **Argonne National Laboratory**;
                   2020;
                   https://www.mcs.anl.gov/petsc
