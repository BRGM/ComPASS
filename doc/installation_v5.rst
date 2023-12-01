.. meta::
    :scope: version5

Installation instructions
=========================

Following `The Zen of Python <https://www.python.org/dev/peps/pep-0020/>`_:

    There should be one -- and preferably only one -- obvious way to do it

Unfortunately subtle differences may exist from one system to another.
Please do not hesitate to submit issues/comments to improve this section.

.. _using conda environments:

Using conda environments
------------------------

As of today, using
`conda environments <https://docs.conda.io/projects/conda/en/latest/user-guide/concepts/environments.html>`_
is probably the fastest way to get ready to use ComPASS.
It will download a lot of dependencies: these may represent a subsequent payload
and will occupy a certain amount of space on your drive. Yet, you will
end up with an isolated environment to use and/or develop ComPASS.
Though ComPASS will still be compiled on your system,
if you are interested in maximising performance and/or exploiting
libraries that have been fine-tuned for your system (e.g. PETSc)
you may want to consider a native build (cf. the next section).

Prerequisites
^^^^^^^^^^^^^

Install conda
"""""""""""""

You will need to have
`conda <https://docs.conda.io/projects/conda/en/latest/index.html>`_
on your system.
Please consider reading the
`installation instructions
<https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html>`_.
The easiest (and lightest) way is probably to install
`Miniconda <https://docs.conda.io/projects/conda/en/latest/glossary.html#miniconda-glossary>`_.
Miniconda installers can be found on this
`download page <https://docs.conda.io/en/latest/miniconda.html#latest-miniconda-installer-links>`_.

Once you have installed conda (either through Miniconda or Anaconda) you should have
the conda command available on your path. You may consider updating conda to have the
latest release.

.. code-block:: bash

  conda update -y conda


Replace conda with mamba (optional)
"""""""""""""""""""""""""""""""""""

`Mamba <https://github.com/mamba-org/mamba>`_
is an implementation of the conda package manager in C++
which is much faster than conda and will speed-up
many of the installation steps.

You can install it from `conda-forge <https://conda-forge.org/>`_.

.. code-block:: bash

  conda activate base
  conda install -c conda-forge mamba

or for one-liners:

.. code-block:: bash

  conda install -n base -c conda-forge mamba

Then, in the following you can replace all the occurences of
the :code:`conda` command with the :code:`mamba` command.


You're almost there
^^^^^^^^^^^^^^^^^^^

In the following, you will need to have conda (or mamba) on your path
and this may require to activate your conda :code:`base` environment
(:code:`conda activate base` or :code:`source activate base` depending
on your settings).

Step 1
""""""

Create the *compass* `conda environment
<https://gitlab.com/compass/compass-v5/compass-sdk/-/blob/f58d4662a4f15a3bedc7bd430570d52d76f95fb2/environment.yml>`_.

.. code-block:: bash

  conda env create -f https://gitlab.com/compass/compass-v5/compass-sdk/-/raw/main/environment.yml

The whole process may be a bit long because many packages will be
downloaded.
Once it is finished, you can activate the *compass* environment with:

.. code-block:: bash

  conda activate compass

In the following we suppose that this *compass* environment has been activated.

Step 2
""""""

Install the latest version of ComPASS.

You should compile and install a `compass utility <https://gitlab.com/compass/compass-v5/compass>`_
using a single line pip command:

.. code-block:: bash

  pip install git+https://gitlab.com/compass/compass-v5/compass.git

Then you can use the *compass* command to install, use, develop
the different software bricks that constitute the ComPASS family of tools.
To list the utilities, run :code:`compass -h`
and/or :code:`compass <command> -h`.

To clone and install all the modules in the *modules* directory,
you need also to define the directory as *root*:

.. code-block:: bash

  mkdir -p modules
  compass set-root ./modules
  compass clone --all
  compass build --all --build-type Release

It is possible to run all the c++ and python tests using
:code:`compass cpptest --all` and :code:`compass pytest --all`.
It might be long as it runs all the tests of all the software bricks.

Once it is done, you can start :ref:`using simulation scripts <Setting-up a simulation>`.

Just **remember to activate the** *compass* **conda environment**
each time you want to use compass.


Additional remarks
^^^^^^^^^^^^^^^^^^

Verbose mode for compass build
""""""""""""""""""""""""""""""

If you need details about the installation step, you can run *compass* in
verbose mode using the :code:`-vvv` flag, just running:

.. code-block:: bash

  compass build --all --build-type Release --verbose


Useful `conda env create` options
""""""""""""""""""""""""""""""""""

:code:`--force`: force creation of environment (removing a previously existing
environment of the same name).

:code:`-n new_name`: will rename the generated environment.

Gitlab servers used for collaborative development
"""""""""""""""""""""""""""""""""""""""""""""""""

*ComPASS* is developed on the
`BRGM gitlab server <https://gitlab.brgm.fr/brgm/modelisation-geologique/compass>`_.
The main branches are mirrored to `gitlab.com <https://gitlab.com/compass/compass-v5>` so that
any gitlab URL above can be replaced with `gitlab.brgm.fr/brgm/modelisation-geologique/`
to have access to development branches.
An access to `BRGM gitlab server <https://gitlab.brgm.fr/brgm/modelisation-geologique/compass>`_
can be provided upon `request <mailto:compass@brgm.fr>`_.

It is possible to clone all the bricks from the BRGM gitlab platform using:

.. code-block:: bash

  compass clone --repo https://gitlab.brgm.fr/brgm/modelisation-geologique/compass --all
