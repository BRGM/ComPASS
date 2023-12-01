Profiling
=========

You can easily profile sequential test using the `built-in python profiler
<https://docs.python.org/3/library/profile.html>`_

Supposing that you have a script that has a timeloop relying
on the `Standard_time_loop` function from the `compass-coats` module.

.. ifconfig:: versionlevel <= '4'

    .. code-block:: python

        import ComPASS
        from ComPASS.timeloops import standard_loop
        simulation = ComPASS.load_physics('water2ph')
        # init your simulation here...
        standard_loop(simulation, some parameters here...)

.. ifconfig:: versionlevel > '4'

    .. code-block:: python

        from compass_coats.models import Coats
        from compass_coats.standard_time_loop import Standard_time_loop
        model = Coats("diphasic")
        # init your simulation here...
        time_loop = Standard_time_loop(
            geom=geom,
            model=model,
            scheme=scheme,
            data=data,
        )
        time_loop.run(
            initial_step=dt,
            final_time=final_time,
        )

You can easily profile only the :code:`time_loop.run` wrapping it
with a few lines of codes:

.. ifconfig:: versionlevel <= '4'

    .. code-block:: python

        import ComPASS
        from ComPASS.timeloops import standard_loop

        simulation = ComPASS.load_physics('water2ph')

        # init your simulation here...

        import cProfile
        import ComPASS.mpi as mpi  # to access mpi.proc_rank

        pr = cProfile.Profile()
        pr.enable()

        standard_loop(simulation, some parameters here...)

        pr.disable()
        pr.dump_stats(f"timeloop-proc{mpi.proc_rank:05d}.profile")

.. ifconfig:: versionlevel > '4'

    .. code-block:: python

        from compass_coats.models import Coats
        from compass_coats.standard_time_loop import Standard_time_loop
        model = Coats("diphasic")
        # init your simulation here...
        time_loop = Standard_time_loop(
            geom=geom,
            model=model,
            scheme=scheme,
            data=data,
        )

        import cProfile
        from compass_utils import mpi  # to access mpi.proc_rank

        pr = cProfile.Profile()
        pr.enable()

        time_loop.run(
            initial_step=dt,
            final_time=final_time,
        )

        pr.disable()
        pr.dump_stats(f"timeloop-proc{mpi.proc_rank:05d}.profile")


Then you can use graphics tools to explore profiling results, such as:
    * `RunSnakeRun <https://pypi.org/project/RunSnakeRun/>`_:
        - `pip install runsnakerun`
        - `runsnake timeloop-proc00000.profile`
    * `SnakeViz <https://jiffyclub.github.io/snakeviz/>`_
