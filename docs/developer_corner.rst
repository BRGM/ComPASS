==================
Developer's corner
==================

Coding conventions
==================

Checkout the CodingConventions.rst file at ComPASS root directory.


Code formatting
===============

From ComPASS 4.2 code formatting is enforced using
the `pre-commit utility <https://pre-commit.com/>`_.
Pipelines will fail for badly formatted code.

Installation
------------

* The ``pre-commit`` and all the required hooks are shipped
  with the ComPASS docker work environment.

  Orherwise, you will need to install ``pre-commit`` and ``clang-format``,
  for example if you are running ubuntu:

  .. code:: shell

        python -m pip install pre-commit
        sudo apt-get install clang-format

  Other dependencies will be installed through ``pre-commit`` installation.

* If you start a new branch from ``v4.2`` or any newer version
  you can skip this point.

  Otherwise go to ComPASS root directory and retrieve the ``pre-commit``
  and ``clang-format`` configuration files.

  .. code:: shell

    git checkout v4.2 -- .pre-commit-config.yaml .clang-format

* Then, install ``pre-commit`` running in ComPASS root directory:

  .. code:: shell

    pre-commit install

You're done...

Using pre-commit
----------------

Once installed, each time you will issue a ``git commit``, ``pre-commit`` will run several hooks
on modified files. The commit will suceed if and only if all files are correctly formatted.

Nevertheless, the hooks will modify the files so that they match format requirements. Eventually you will
just have to add the modifications (with ``git add`` or the more risky ``git commit -a``)
and try to issue the commit again. It should succeed.

If you want to manually run ``pre-commit`` hooks on a specifc set of files you can use the
run subcommand (cf. the `command documentation <https://pre-commit.com/#pre-commit-run>`_):

    .. code:: shell

        pre-commit run --file path1 path2


.. note::
    Most of formatting choices were set to the hooks default.
    Don't hesitate to create new gitlab issues to discuss formatting choices.


Reformat branch history
-----------------------

A script to duplicate a branch applying formating tools
is available in ``sdk/format_history.py``.

.. warning::
    For the time being you will need python 3.7 or higher.
    All commits on the ``feature`` branch must contain code that is syntaxically corect.
    Otherwise the formatting hooks may fail.
    
To display the help message, run:

.. code:: shell

    python sdk/format_history.py --help


The common use is to reformat a branch (let's say ``feature``)
that has diverged from a previous version (let's say ``v4.1``)
before submitting a merge request to a newer version (let's say ``v4.2``).

Then running:

.. code:: shell

    python sdk/format_history.py v4.1 feature --format-source v4.2

will first issue a commit with the *"Reformat code"* message that will contain
a formatted version of all the files present in `v4.1` that were modified in `feature`.

Then all the commits in the ``feature`` branch will be re-issued with a formatted version of the files
in a *detached* HEAD state.

Then if the reformatting script works fine, you can save the result in a new branch, *e.g.*:

.. code:: shell

    git checkout -b reformatted_feature


Finally you are ready to rebase ``reformatted_feature`` on ``v4.2`` skipping the
first *"Reformat code"* commit. Supposing that this commit has sha1 ``abcdef0``,
you could try something like :

.. code:: shell

    git rebase abcdef0 reformatted_feature --onto v4.2


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
