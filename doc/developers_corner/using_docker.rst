.. _using_compass_sdk_with_docker:

Using *ComPASS SDK* via docker images
=====================================

Though these docker environments are still used for Continuous Integration
their use is discouraged in favor of conda environments (cf. :ref:`using conda environments`).

Selecting the appropriate docker image
--------------------------------------

The version of the *ComPASS SDK* that is used in CI Pipelines
can be found at the top of the `.gitlab-ci.yml` file
located in the root directory.
It is the value of the `COMPASS_SDK` variable.

The docker images are generated via a
`separate gitlab project <https://gitlab.inria.fr/charms/compass-images>`_
and can be found under the `Container Registry header <https://gitlab.inria.fr/charms/compass-images/container_registry>`_.
The value of the `COMPASS_SDK` variable can then be used
to select the appropriate docker image.

Work environment
----------------

The `charms/compass-mages/work <https://gitlab.inria.fr/charms/compass-images/container_registry/1149>`_
containers are designed to
**compile and modify ComPASS and run simulations**.
They ship all what is needed: compilation environment
with dependencies and a few goodies such as scipy or matplotlib modules.
Consequently they are relatively big containers but
you will only need to install (one of) them once!
The idea is to work on your local/host file system through the ``/localfs``
container volume (\ ``localfs``\ for *local file system*\ ).

Let suppose that you work in a local directory called ``/my/local/dir``.
If not already done clone the ComPASS git repository:

.. code-block:: shell

   cd /my/local/dir
   git clone git@gitlab.inria.fr:charms/ComPASS.git

or via https if you prefer:

.. code-block:: shell

   cd /my/local/dir
   git clone https://gitlab.inria.fr/charms/ComPASS.git

Once you have installed docker you can pull the container
from the registry (log in with you gitlab id and password):

.. code-block:: shell

   docker login registry.gitlab.inria.fr
   docker pull registry.gitlab.inria.fr/charms/compass-images/work:0.1.5

Then you can possibly tag it to use shorter name:

.. code-block:: shell

   docker tag registry.gitlab.inria.fr/charms/compass-images/work:0.1.5 compass-work

Container are created with a default user which is called *compass* and has uid 1000.
When mounting a a volume you may experience problems if your host uid is not the user uid.
In this case use the ``--user`` flag from the docker run command
(cf. `docker documentation <https://docs.docker.com/engine/reference/run/#user>`_\ ).
Now you can run the container interacting with the code that is on your local filesystem
mouting ``/my/local/dir`` through the ``/localfs`` volume:

.. code-block:: shell

   docker run -it -v /my/local/dir:/localfs compass-work --compass-uid $UID

Depending on your platform the variable ``UID`` may not be defined
and you can replace it with ```id -u $USER` ``.

Check that you see everything that you have in ``/my/local/dir`` through the ``/localfs``
of the container, i.e. once in the container ``ls /localfs``
should list the content of ``/my/local/dir`` and esspecially the ``ComPASS`` source directory.

Now, inside the container you can just compile and run ComPASS
(let's suppose that you have 4 cores on you host you can use them on the container):

.. code-block:: shell

   cd /localfs/ComPASS
   python3 setup.py develop --build-type=Debug -G "Unix Makefiles" -j 4 -DComPASS_WITH_ALL_PHYSICS=ON
   cd test/bulk
   mpirun -n 4 python3 doublet_on_cartesian_grid.py

All changes that you make in the source files will remain on your local file system
(and you may want to push them to dedicated branches on the git repository).
