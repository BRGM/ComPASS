==================
Developer's corner
==================

Coding conventions
==================

Checkout the CodingConventions.rst file at ComPASS root directory.


Profiling
=========

You can easily profile sequential test using the `built-in python profiler
<https://docs.python.org/3/library/profile.html>`_

Supposing that you have a script that has a timeloop relying on the `standard_loop` 
function from the `timeloops` module.

.. code-block:: python

    import ComPASS
    from ComPASS.timeloops import standard_loop

    simulation = ComPASS.load_eos('water2ph')

    # init your simulation here...

    standard_loop(simulation, some parameters here...)

You can easily profile only the standard_loop wrapping it with a few lines of codes:

.. code-block:: python

    import ComPASS
    from ComPASS.timeloops import standard_loop
    from ComPASS import mpi # to access mpi.proc_rank

    simulation = ComPASS.load_eos('water2ph')

    # init your simulation here...

    import cProfile
    import ComPASS.mpi as mpi

    pr = cProfile.Profile()
    pr.enable()

    standard_loop(simulation, some parameters here...)

    pr.disable()
    pr.dump_stats(f"timeloop-proc{mpi.proc_rank:05d}.profile")


Then you can use graphics tools to explore profiling results, such as:
    * `RunSnakeRun <https://pypi.org/project/RunSnakeRun/>`_:
        - `pip install runsnakerun`
        - `runsnake timeloop-proc00000.profile`
    * `SnakeViz <https://jiffyclub.github.io/snakeviz/>`_
