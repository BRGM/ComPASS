.. meta::
    :scope: version5

Outputs
=======

The output files are in the directory specified when creating the instance
of the :code:`compass-coats.Standard_time_loop`.
By default the directory is named *output*.
It is recommended to specify a unique name
to avoid erasing accidently some results.
Using :code:`output_directory(__file__)`
stores the output files in the *output-{file_name}* directory.

.. code-block:: python

    solution, tick = time_loop.run(
        initial_step=dt,
        final_time=final_time,
        output_period=snapshot_period,
        timeloop_file=True,
    )

The outputs contain :

* states at particular time steps (those specified with *output_period* and/or
  *output_every* and/or *nb_output*). To prepare them for visualization, execute

  .. code-block:: python

    postprocess(output_dir, time_unit="day")

  By default output_dir="output" and time_unit="year". *time_unit* allows
  second, minute, hour, day or year.

* information about the convergence behaviour (if *timeloop_file=True*).
  The yaml file *timeloop_log.yaml* contains a dictionary with all
  the time steps,
  for each time step it contains convergence information such as
  the computation time,
  the number of Newton iterations, the number of linear solver iterations
  for each Newton iteration.
