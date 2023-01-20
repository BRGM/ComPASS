.. _development_projects:

Credits
=======

ComPASS is co-developped by `BRGM`_ and
`Université Côte d'Azur`_ (`LJAD`_ - `Inria`_)
and the main evolutions were/are implemented through the following projects:

.. admonition:: 2021 - ...... |sp| Partnership with `ANDRA`_ |brsp| Application to radioactive waste storage
  :class: curriculum-vitae

  * Cell-centered discretization
  * Improve modularity of the code
  * Prospective research on Thermo Hydro Mechanical models

.. admonition:: 2019 - 2022 |sp| Collaboration with `Storengy`_ |brsp| Thermal well modelling for geothermal systems
  :class: curriculum-vitae

  * Multi-branch wells
  * Drift Flux Model
  * Geothermal test cases

.. admonition:: 2017 - 2020 |sp| `ANR CHARMS`_ with `BRGM`_, `Storengy`_, `MdlS`_, `LJLL`_, `Université Côte d'Azur`_
  :class: curriculum-vitae

  * High level Python API
  * Atmospheric boundary condition
  * Benchmark and geothermal tests
  * HPC profiling with Scalasca and Paraver (workshop EOCOE)

  .. admonition:: L. Beaude, K. Brenner, S. Lopez, R. Masson, F. Smai,
    *Non-isothermal compositional liquid gas darcy flow: formulation,
    soil-atmosphere boundary condition and application to high energy
    geothermal simulations*,
    `Computational Geosciences, 2019 <https://doi.org/10.1007/s10596-018-9794-9>`_
    :class: publication

    | Computational Geosciences
    | Volume 23,
    | 2019,
    | Pages 443-470,
    | **Abstract:** This article deals with the modelling and formulation of
      compositional gas liquid Darcy flow. Our model includes an advanced
      boundary condition at the interface between the porous medium and the
      atmosphere accounting for convective mass and energy transfer, liquid
      evaporation and liquid outflow. The formulation is based on a fixed set
      of unknowns whatever the set of present phases. The thermodynamic
      equilibrium is expressed as complementarity constraints. The model and
      its formulation are applied to the simulation of the Bouillante
      high-energy geothermal field in Guadeloupe characterised by a high
      temperature close to the surface.
    | **Keywords:** Non-isothermal compositional Darcy flow; Geothermal energy;
      Soil-atmosphere boundary condition; Outflow boundary condition;
      Porous medium drying; Finite volume scheme

.. admonition:: 2016 |sp| CEMRACS summer camp (`BRGM`_, `Université Côte d'Azur`_, `CEA`_)
  :class: curriculum-vitae

  * First well model implementation

  .. admonition:: L. Beaude, T. Beltzung, K. Brenner, S. Lopez, R. Masson,
    F. Smai, J.F. Thebault, F. Xing,
    *Parallel geothermal numerical model with fractures and multi-branch wells*,
    `ESAIM: Proceedings, 2018 <https://doi.org/10.1051/proc/201863109>`_
    :class: publication

    | ESAIM: Proceedings and surveys,
    | Volume 63,
    | 2018,
    | CEMRACS 2016,
    | Pages 109-134,
    | **Abstract:** To answer the need for an efficient and robust geothermal
      simulation tool going beyond existing code capabilities in terms of
      geological and physical complexity, we have started to develop a
      parallel geothermal simulator based on unstructured meshes. The model
      takes into account complex geology including fault and fracture
      networks acting as major heat and mass transfer corridors and complex
      physics coupling the mass and energy conservations to the thermodynamic
      equilibrium between the gas and liquid phases. The objective of this
      Cemracs project was to focus on well modeling which is a key missing
      ingredient in our current simulator in order to perform realistic
      geothermal studies both in terms of monitoring and in terms of history
      matching. The well is discretized by a set of edges of the mesh in
      order to represent efficiently slanted or multi-branch wells on
      unstructured meshes. The connection with the 3D matrix and the 2D
      fracture network at each node of the well is accounted for using
      Peaceman's approach. The non-isothermal flow model inside the well is
      based on the usual single unknown approach assuming the hydrostatic and
      thermodynamical equilibrium inside the well. The parallelization of the
      well model is implemented in such a way that the assembly of the
      Jacobian at each Newton step and the computation of the pressure drops
      inside the well can be done locally on each process without MPI
      communications.


.. admonition:: 2015 - 2016 |sp| inter-Carnot `BRGM`_-`Inria`_ postdoctoral fellowship of Feng Xing
  :class: curriculum-vitae

  * Multi-phase multi-component non-isothermal flow in fractured porous media for high energy geothermal systems

  .. admonition:: F. Xing, R. Masson, S. Lopez, *Parallel numerical modeling of
    hybrid-dimensional compositional non-isothermal Darcy flows in fractured porous media*,
    `JCP, 2017 <https://doi.org/10.1016/j.jcp.2017.05.043>`_
    :class: publication

    | Journal of Computational Physics
    | Volume 345,
    | 2017,
    | Pages 637-664,
    | ISSN 0021-9991,
    | **Abstract:** This paper introduces a new discrete fracture model accounting for non-isothermal compositional
      multiphase Darcy flows and complex networks of fractures with intersecting, immersed and non-immersed fractures.
      The so called hybrid-dimensional model using a 2D model in the fractures coupled with a 3D model in the matrix
      is first derived rigorously starting from the equi-dimensional matrix fracture model. Then, it is discretized
      using a fully implicit time integration combined with the Vertex Approximate Gradient (VAG) finite volume scheme
      which is adapted to polyhedral meshes and anisotropic heterogeneous media. The fully coupled systems are assembled
      and solved in parallel using the Single Program Multiple Data (SPMD) paradigm with one layer of ghost cells.
      This strategy allows for a local assembly of the discrete systems. An efficient preconditioner is implemented to
      solve the linear systems at each time step and each Newton type iteration of the simulation. The numerical
      efficiency of our approach is assessed on different meshes, fracture networks, and physical settings in terms of
      parallel scalability, nonlinear convergence and linear convergence.
    | **Keywords:** Discrete fracture model; Non-isothermal compositional multiphase hybrid-dimensional Darcy flow
      model; Vertex Approximate Gradient scheme; Polyhedral meshes; Parallel algorithm; Preconditioner

  .. admonition:: F. Xing, R. Masson, S. Lopez, *Parallel Vertex Approximate
    Gradient discretization of hybrid dimensional Darcy flow and transport
    in discrete fracture networks*,
    `Computational Geosciences, 2016 <https://doi.org/10.1007/s10596-016-9606-z>`_
    :class: publication

    | Computational Geosciences
    | Volume 21,
    | 2016,
    | Pages 595-617,
    | **Abstract:** This paper proposes a parallel numerical algorithm to
      simulate the flow and the transport in a discrete fracture network taking
      into account the mass exchanges with the surrounding matrix. The
      discretization of the Darcy fluxes is based on the Vertex Approximate
      Gradient finite volume scheme adapted to polyhedral meshes and to
      heterogeneous anisotropic media, and the transport equation is
      discretized by a first-order upwind scheme combined with an Euler
      explicit integration in time. The parallelization is based on the single
      program, multiple data (SPMD) paradigm and relies on a distribution of
      the mesh on the processes with one layer of ghost cells in order to allow
      for a local assembly of the discrete systems. The linear system for the
      Darcy flow is solved using different linear solvers and preconditioners
      implemented in the PETSc and Trilinos libraries. The convergence of the
      scheme is validated on two original analytical solutions with one and
      four intersecting fractures. Then, the parallel efficiency of the
      algorithm is assessed on up to 512 processes with different types of
      meshes, different matrix fracture permeability ratios, and different
      levels of complexity of the fracture network.


.. admonition:: 2012 - 2014 |sp| C. Guichard, R. Eymard et al. (`Université Côte d'Azur`_, `IFPEN`_)
  :class: curriculum-vitae

  * Parallel implementation of the Vertex Approximate Gradient scheme
  * First two-phase Darcy flow implementation

  .. admonition:: R. Eymard, C. Guichard, R. Masson,
    *High Performance Computing linear algorithms for two-phase flow
    in porous media*,
    `FVCA7, 2014 <https://doi.org/10.1007/978-3-319-05591-6_55>`_
    :class: publication

    | Finite Volumes for Complex Applications VII-Elliptic, Parabolic and
      Hyperbolic Problems. Springer Proceedings in Mathematics & Statistics,
    | Volume 78,
    | 2014,
    | Pages 557-565,
    | **Abstract:** We focus here on the difficult problem of linear solving,
      when considering implicit scheme for two-phase flow simulation in
      porous media. Indeed, this scheme leads to ill-conditioned linear
      systems, due to the different behaviors of the pressure unknown (which
      follows a diffusion equation) and the saturation unknown (mainly
      advected by the total volumic flow). This difficulty is enhanced by the
      parallel computing techniques, which reduce the choice of the possible
      preconditioners. We first present the framework of this study, and then
      we discuss different algorithms for linear solving. Finally, numerical
      results show the performances of these algorithms.

  .. admonition:: E. Dalissier, C. Guichard, P. Havé, R. Masson, C. Yang,
    *ComPASS: a tool for distributed parallel
    finite volume discretizations on general unstructured polyhedral meshes*,
    `ESAIM: Proceedings, 2013 <https://doi.org/10.1051/proc/201343010>`_
    :class: publication

    | ESAIM: Proceedings,
    | Volume 43,
    | 2013,
    | CEMRACS 2012,
    | Pages 147-163,
    | **Abstract:** The objective of the ComPASS project is to develop a
      parallel multiphase Darcy flow simulator adapted to general unstructured
      polyhedral meshes (in a general sense with possibly non planar faces) and
      to the parallelization of advanced finite volume discretizations with
      various choices of the degrees of freedom such as cell centres, vertices,
      or face centres. The main targeted applications are the simulation of CO2
      geological storage, nuclear waste repository and reservoir simulations.
      The CEMRACS 2012 summer school devoted to high performance computing has
      been an ideal framework to start this collaborative project. This paper
      describes what has been achieved during the four weeks of the CEMRACS
      project which has been focusing on the implementation of basic features
      of the code such as the distributed unstructured polyhedral mesh, the
      synchronization of the degrees of freedom, and the connection to
      scientific libraries including the partitioner METIS, the visualization
      tool PARAVIEW, and the parallel linear solver library PETSc. The parallel
      efficiency of this first version of the ComPASS code has been validated
      on a toy parabolic problem using the Vertex Approximate Gradient finite
      volume spatial discretization with both cell and vertex degrees of
      freedom, combined with an Euler implicit time integration.

.. papiers !

.. _ANDRA: https://www.andra.fr/
.. _ANR CHARMS: http://www.anr-charms.org/
.. _BRGM: https://www.brgm.fr/fr
.. _CEA: https://www.cea.fr/
.. _Engie: https://www.engie.fr/
.. _IFPEN: https://www.ifpenergiesnouvelles.fr/
.. _Inria: https://www.inria.fr/fr/centre-inria-universite-cote-azur
.. _LJAD: https://univ-cotedazur.fr/laboratoires/laboratoire-jean-alexandre-dieudonne-ljad
.. _LJLL: https://www.ljll.math.upmc.fr/
.. _MdlS: https://mdls.fr/
.. _Storengy: https://www.storengy.com/fr
.. _Université Côte d'Azur: https://univ-cotedazur.fr/

.. |brsp| raw:: html

  <br/>&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp

.. |sp| raw:: html

  &nbsp&nbsp&nbsp
