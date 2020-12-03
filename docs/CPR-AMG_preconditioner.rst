The CPR-AMG preconditioner in ComPASS
=======================================

Introduction
-----------------

Linear systems in ComPASS describe complex physics on mesh degrees of
freedom of different natures (pressure, temperature and wells), governed
by partial differential equations with distinct analytical properties.
Usual preconditioning methods like ILU, Schwarz methods and multigrid
algorithms prove to be unsufficient on large systems because they treat
all unknowns the same way, thus losing knowledge of their physical and
numerical properties. The default preconditioner in ComPASS is CPR-AMG,
a method which uses the splitting of the blocks of unknowns in the linear
operator to improve the convergence of Krylov methods. This document
describes the CPR-AMG preconditioning procedure and its Python
implementation in the framework of the linear algebra package for
ComPASS.

The CPR-AMG preconditioning method
----------------------------------

Let :math:`A` be a linear operator describing the behaviour of type of
unknowns : pressure and temperature.

.. math::

   A=\left(\begin{array}{ll}
   A_{p p} & A_{pT} \\
   A_{T p} & A_{T T}
   \end{array}\right)

We are looking to solve the following linear system : :math:`A x=b`, with
:math:`x=\left(\begin{array}{l} x_{p} \\ x_{T} \end{array}\right) \quad
\text{and } b=\left(\begin{array}{l} b_{p} \\ b_{T} \end{array}\right)`,
using the CPR (Constrained Pressure Residual) preconditioner
([Roy_2019]_). Application of this preconditioner to a
vector :math:`b` is performed like so :

.. math::

      \left\{\begin{array}{l}
      y=\left(\begin{array}{l}
      y_{p} \\
      y_{T}
      \end{array}\right)=\left(\begin{array}{l}
      0 \\
      0
      \end{array}\right) \\
      \tilde{A}_{p p} y_{p}=b_{p} \\
      r=b-A y \\
      \tilde{L} \tilde{U} x=r \\
      x=x+y
      \end{array}\right.

First a partial resolution is performed on the pressure block, then a
correction is made to the original vector :math:`b`. Finally a
pseudo-solve is applied to the full set of unknowns. In the CPR-AMG
method, :math:`A_{p p}^{-1}` is approximated by an algebraic multigrid
v-cycle, and :math:`A^{-1}` by an incomplete LU solve, denoted by the
:math:`\sim` markers. This procedure corresponds to the multiplication
by this matrix :

.. math::

   \text{CPR}(A) = (\tilde{L}\tilde{U})^{-1}\left(I-A M\right)+M

with
:math:`M=\left(\begin{array}{cc} \tilde{A}_{p p}^{-1} & 0 \\ 0 & 0 \end{array}\right)`.

Implementation with petsc4py
----------------------------

Preconditioners of the form
:math:`M^{-1} = M_2^{-1}\left(I-A M_1^{-1}\right)  + M_1^{-1}` like CPR
are implemented under the type ``COMPOSITE``
in multiplicative mode (additive application corresponds to the matrix
:math:`M^{-1} = M_2^{-1} + M_1^{-1}`). Documentation can be found at
page 90 in [PETSc_manual]_.

.. literalinclude:: ../ComPASS/linalg/petsc_linear_solver.py
   :language: python
   :pyobject: PetscIterativeSolver.set_cpramg_pc
   :start-at: cpramg_pc = self.pc
   :end-before: setUp

| The ``PC`` objects must be instanciated in the required order. In the
  case of CPR-AMG, the ``PC`` :math:`M_2` is an ILU resolution (default
  in PETSc), defined by diagonal blocks with the ``BJACOBI`` option for
  parallel compatibility. The ``PC`` :math:`M_1` is set as a
  ``FIELDSPLIT`` type preconditioner ([PETSc_manual]_
  page 93). This type is a special case of ``COMPOSITE`` where sub-PCs
  are restricted to selected sets of unknowns called "fields". Here the
  fields are set using the IS (Index Sets) objects in PETSc. These
  consist in integer arrays which map the indices of the local
  pressure submatrix with the corresponding global indices of the
  pressure unknowns in the global matrix. (e.g. Let :math:`N_p` be the number
  of pressure unknowns on current processor :math:`p`. Then for :math:`i \in
  [0, Np-1]`, :math:`\text{p_IS}[i]` is the global index of the
  local :math:`i`-th pressure unknown in matrix :math:`A`). Note that the
  well unknowns are included in the pressure field, because we want them
  to be affected by the AMG procedure together with the pressure unknowns.
  Then the IndexSet for the remaining unknowns is built using the PETSc
  function ``complement``, (e.g. every index which is not in the pressure
  IS).
| The type of the fieldsplit ``PC`` is set on ``ADDITIVE``, e.g.
  :math:`M_{FS} = M_{FSp} + M_{FST}` with :math:`M_{FST} = 0`. This type
  of preconditioner is not defined in PETSc, meaning we have to set it
  ourselves with a custom PC class which returns zero.

.. literalinclude:: ../ComPASS/linalg/petsc_linear_solver.py
   :language: python
   :pyobject: PetscIterativeSolver.set_cpramg_pc
   :start-at: # The Index Set for the pressure field
   :end-at: null_pc.setPythonContext(NullPC())

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
