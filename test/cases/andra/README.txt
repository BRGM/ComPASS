 * Cigeo simulation (Mira Kabbara):
some factories for the mesh are in the maillage_Cigeo directory.
It contains the file advanced-extract-mesh.py which extracts a sub-mesh
containing one gallery, removes the shrinked cells and identifies
the necessary nodes. You can also visualize the meshes using the vtu format.

First isothermal simulation on the small mesh:
    - init step with no cell in the gallery and atmospheric BC at the gallery
    walls to init with some gas inside
    python small_cigeo_isotherm_atm_gal_init.py
    - then reload and execute Cigeo simulation
    python small_cigeo_isotherm_reload.py

An anisotropic bulk thermal conductivity test is available with:
     - the initial stage (atmospheric BC instead of the gallery)
     python small_cigeo_anisotropic_atm_gal_init.py
     - the Cigeo simulation after reloading
     python small_cigeo_anisotropic_reload.py

The same anisotropic test with Neumann heat flux representing the waste
is under preparation with:
     - the same initial stage (atmospheric BC instead of the gallery)
     python small_cigeo_anisotropic_atm_gal_init.py
     - the Cigeo simulation after reloading
     python small_cigeo_anisotropic_reload_Neumann_heat_flux.py


 * Nabil's test cases:
one gallery in a 3D non-cartesian mesh,
the mesh is small_disposal_gallery.msh
python andra_small_disposal_gallery_gmsh.py
