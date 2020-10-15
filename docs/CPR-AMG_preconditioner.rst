The CPR-AMG preconditioner in ComPASS
================================

Introduction
-----------------

Linear systems in ComPASS describe complex physics on mesh degrees of
freedom of different natures (pressure, temperature and wells), governed
by partial differential equations with distinct analytical properties.
Usual preconditioning methods like ILU, Schwarz methods and multigrid
algorithms prove to be unsufficient on large systems because they treat
all unknowns the same way, thus losing knowledge of their physical and
numerical properties. The default preconditioner in ComPASS is CPR-AMG,
a method which uses the knowledge of the unknowns blocks in the linear
operator to improve the convergence of Krylov methods. This document
describes the CPR-AMG preconditioning procedure and its Python
implementation in the framework of the linear algebra package for
ComPASS.

The CPR-AMG preconditioning method
-----------------

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
   \tilde{A}_{p p} x_{p}=b_{p} \\
   b'_{T}=b_{T}-A_{T p} x_{p} \\
   \tilde{A} y=b'
   \end{array}\right.

First a partial resolution is performed on the pressure block, then a
correction is made to the original vector :math:`b`. Finally a
pseudo-solve is applied to the full set of unknowns. In the CPR-AMG
method, :math:`A_{p p}^{-1}` is approximated by an algebraic multigrid
v-cycle, and :math:`A^{-1}` by an incomplete LU solve, denoted by the
:math:`\sim` markers. This procedure corresponds to the multiplication
by this matrix :

.. math::

   \text{CPR}(A) = \tilde{A}^{-1}\left(I-A M\right)+M

with
:math:`M=\left(\begin{array}{cc} \tilde{A}_{p p}^{-1} & 0 \\ 0 & 0 \end{array}\right)`.

Implementation with petsc4py
-----------------

Preconditioners of the form
:math:`M^{-1} = M_2^{-1}\left(I-A M_1^{-1}\right)  + M_1^{-1}` like CPR
are implemented under the type ``COMPOSITE``
in multiplicative mode (additive application corresponds to the matrix
:math:`M^{-1} = M_2^{-1} + M_1^{-1}`). Documentation can be found at
page 90 in [PETSc_manual]_.

::

     import petsc4py
     petsc4py.init(sys.argv)
     from petsc4py import PETSc

     cpramg_pc = PETSc.PC().create(comm=comm)
     cpramg_pc.setOperators(A, A)
     cpramg_pc.setType(PETSc.PC.Type.COMPOSITE)
     cpramg_pc.setCompositeType(PETSc.PC.CompositeType.MULTIPLICATIVE)
     cpramg_pc.addCompositePC(PETSc.PC.Type.FIELDSPLIT)
     cpramg_pc.addCompositePC(PETSc.PC.Type.BJACOBI)

| The ``PC`` objects must be instanciated in the required order. In the
  case of CPR-AMG, the ``PC`` :math:`M_2` is an ILU resolution (default
  in PETSc), defined by diagonal blocks with the ``BJACOBI`` option for
  parallel compatibility. The ``PC`` :math:`M_1` is set as a
  ``FIELDSPLIT`` type preconditioner ([PETSc_manual]_
  page 93). This type is a special case of ``COMPOSITE`` where sub-PCs
  are restricted to selected sets of unknowns : "``fields``". It is
  possible to distinguish these sets specifically using the ``IS``
  (index sets) object type in PETSc, but if the unknowns are stored in
  blocks, which is the case of the matrices from ComPASS, PETSc will
  instanciate these objects internally from the provision of the blocks
  size (physics dependent in ComPASS) and the position of the different
  unknowns in each block (pressure is always the first of each block in
  ComPASS).
| The type of the fieldsplit ``PC`` is ``ADDITIVE``, e.g.
  :math:`M_{FS} = M_{FSp} + M_{FST}` with :math:`M_{FST} = 0`, set with
  type ``PETSc.PC.Type.NONE``

::

     # Split blocks into pressure and the rest of unknowns
     p_field = np.array([0], dtype='int32')
     rest = np.arange(1, block_size, dtype='int32')
     fs_pc = cpramg_pc.getCompositePC(0)
     fs_pc.setFieldSplitFields(block_size, ('pressure', p_field), ('rest of unknowns', rest))
     fs_pc.setFieldSplitType(PETSc.PC.CompositeType.ADDITIVE)
     sub_ksp_list = fs_pc.getFieldSplitSubKSP() # Default KSP is PREONLY in a FieldSplit PC

     # Algebraic multigrid on the pressure field
     pressure_ksp = sub_ksp_list[0]
     pressure_pc = pressure_ksp.getPC()
     pressure_pc.setType(PETSc.PC.Type.HYPRE)
     PETSc.Options().setValue("-pc_hypre_boomeramg_strong_threshold", "0.5")


     # NONE PC on the rest of the unknowns
     rest_ksp = sub_ksp_list[1]
     rest_ksp.getPC().setType(PETSc.PC.Type.NONE)

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
