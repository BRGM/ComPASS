Install troubleshooting
=======================

OpenMPI 4 and mpi4py
--------------------

This issue is only relevant if you use OpenMPI-4 **and** a release of mpi4py
older than 3.0.2. Unfortunately, PETSc 3.13 still installs mpi4py 3.0.1.

OpenMPI 4 now forbids the use of some old MPI-1 features that were removed in MPI-3
(MPI-UB and MPI-LB). This has been fixed in release 3.0.2. Otherwise a workaround
is to modify the file ``src/lib-mpi/config/openmpi.h`` if (installed with PETSc,
it is found in ``$(PETSC_DIR)/$(PETSC_ARCH)/externalpackages/mpi4py-3.0.1``)

A patch is as follows:

.. code:: shell

    *** 142,153 ****
      #endif

      #if OMPI_NUMVERSION >= 40000
    - #ifndef MPI_LB
      #undef  PyMPI_HAVE_MPI_LB
    - #endif
    - #ifndef MPI_UB
      #undef  PyMPI_HAVE_MPI_UB
    - #endif
      #endif /* OMPI >= 4.0.0 */

      #endif /* !PyMPI_CONFIG_OPENMPI_H */
    --- 142,149 ----
