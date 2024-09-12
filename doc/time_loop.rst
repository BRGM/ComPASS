.. meta::
    :scope: version5

Time loop execution
===================

The time loop is usualy an instance of the
:code:`compass-coats.Standard_time_loop` class.

The minimal use is to specify:

* *geom* with the geometry utilities and the mesh,

* *model* with the model and the physics,

* *scheme* with the numerical scheme definition,

* *data* with the petrophysics, the boundary conditions, the initial states...

It is convenient to also specify:

* *output_dir* with the output directory (where
  the visualization files are written), by default it is named *output*,

* *verbosity* (between 0 and 3, where 0 is the less verbose), 3 by default

.. code-block:: python

   from compass_coats.standard_time_loop import Standard_time_loop
   from compass_coats.output_visu import output_directory
   visu_dir = output_directory(__file__)  # is output-{file_name}
   time_loop = Standard_time_loop(
      geom=geom,
      model=model,
      scheme=scheme_def,
      data=data,
      output_dir=visu_dir,
      verbosity=3,
   )

Then to run the time loop, use the :code:`time_loop.run` method.
It allows lots of optional values to set output options or
time-step specifications or callbacks...

.. code-block:: python

   solution, tick = time_loop.run(
      initial_step=0.1 * year,
      final_time=10 * year,
      output_period=final_time / 20,
      output_every=10,  # outputs every 10 iterations
   )
