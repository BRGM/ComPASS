=====================================
Using the ComPASS docker environments
=====================================

Quick links:

* :ref:`Installation instructions<installation>`
* :ref:`Run ComPASS simulations<execute>`
* :ref:`Develop ComPASS<work>`

.. _installation:

Install docker on your machine
==============================

On Ubuntu
^^^^^^^^^

3 steps:


* install,
* add yourself to the docker group
* reboot

.. code-block:: bash

   sudo apt-get install docker.io
   sudo gpasswd -a $USER docker
   sudo reboot

You may want to use the latest docker version to use build steps aliases
(From XXX as XXX), then you can replace the line :

.. code-block:: bash

   sudo apt-get install docker.io

with:

.. code-block:: bash

   curl -sSL https://get.docker.com/ | sh

On Windows
^^^^^^^^^^

go there https://store.docker.com/editions/community/docker-ce-desktop-windows and follow instructions

On Windows 7
~~~~~~~~~~~~

On windows7 docker will run on an underlying virtual machine (let's call it *default*\ ).
This is an important point as connection problems (e.g. internet access via proxies may need a specific configuration direcly on the virtual machine) or mapping between local filesystem and container volumes.

With an already installed version of Oracle VirtualBox accessible on your path, you may want to install a minimum set of additional tools:


#. download the latest version of docker-machine from https://github.com/docker/machine/releases/
#. once docker-machine is on your path, create a *default* virtual machine: the following line will create a virtual machine with 4 cpu and 8GB of RAM

.. code-block:: shell

   docker-machine create --driver virtualbox --virtualbox-cpu-count "4"  --virtualbox-memory "8096" default


#. setup the environment configuration so that docker is accessible: the following command will tell you what to do depending on your shell (here with DOS ``cmd``\ )

.. code-block:: shell

   docker-machine env default --shell cmd

Copy/paste en execute the last line it outputs


#. check your machine is active and running

.. code-block:: shell

   docker-machine ls

The previous command will display the name (here *default*\ ) and the URL that you can use to connect to the machine. You can also use:

.. code-block:: shell

   docker-machine ip machine_name

Where ``machine_name`` is your machine name (here *default*\ )

Then to connect to the machine you can :


* use ``docker-machine ssh machine_name`` (\ ``machine_name`` will default to *default* machine created so you can neglect it if you have created only one virtual machine),
* ssh using the given IP as ``docker`` user (default password ``tcuser``\ ).

You're ready to `use <#execute-compass>`_ docker.

Changing the keyboard layout
""""""""""""""""""""""""""""

When connecting to boot2docker you may end up with a qwerty keyboard layout.

To change it use the solution `here <https://stackoverflow.com/questions/31327923/change-keyboard-layout-boot2docker-tinycore>`_.

When logged as ``docker`` user (default password ``tcuser``\ ) do:

.. code-block:: bash

   tce-load -wi kmaps
   sudo loadkmap < /usr/share/kmap/azerty/fr-latin9.kmap

You may need `this <https://commons.wikimedia.org/wiki/File:QWERTY_keyboard_diagram.svg>`_ to enter the two lines above...

Time synchronization problem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sometimes host and clock lose synchronization and this may lead to spurious messages at buil time.
You can force the time of the VM using:

.. code-block:: shell

   docker-machine ssh default "sudo date -u $(date -u +%m%d%H%M%Y)"

cf also this `link <https://stackoverflow.com/questions/24551592/how-to-make-sure-dockers-time-syncs-with-that-of-the-host>`_.

On MacOS
^^^^^^^^

Follow the instructions available here https://docs.docker.com/docker-for-mac/install/

.. _execute:

Execute ComPASS through docker
==============================

Basic usage
-----------

Login to the gitlab docker registry (cf. `Packages icon on the left <https://gitlab.inria.fr/charms/ComPASS/container_registry>`_\ )

.. code-block:: shell

   docker login registry.gitlab.inria.fr

If you only want to retrieve the latest (develop) container image just pull it:

.. code-block:: shell

   docker pull registry.gitlab.inria.fr/charms/compass:develop

Now the latest version should be among your local docker container images that can be listed with:

.. code-block:: shell

   docker image ls

In the following we suppose that you gave your conainer a shorter name using an alias:

.. code-block:: shell

   docker tag registry.gitlab.inria.fr/charms/compass:develop compass

Then, if you list again all images available (\ ``docker image ls``\ ), you should see the original image and the alias.

The aliased container can now be used as a command whose options can be listed through the ``-h`` or ``--help`` flags:

.. code-block:: shell

   docker run -it compass -h

The idea is to have simulation scripts on your host file system and mount the directory where they lie as a volume when starting the container.
Now you can run the container interacting with scripts that are on your local filesystem mouting ``/my/local/dir`` through the ``/localfs`` volume
with the `-v <https://docs.docker.com/engine/reference/commandline/run/#mount-volume--v---read-only>`_ ``docker run`` mapping option :

.. code-block:: shell

   docker run -it -v /my/local/dir:/localfs compass

The previous command will launch the python interpreter that is used to run script (i.e. from which you can import the ComPASS python module).

The following example suppose that the script ``vertical_column.py`` is in your current working directory
(cf. the ``$PWD`` variable passed to the -v mapping option).

.. code-block:: shell

   docker run -it -v $PWD:/localfs compass vertical_column.py

You may then postprocess outputs using the postprocess script.
Beware that options passed to the postprocess script must be passed using the MS dos way,
i.e. prepending them with a slah (/) instead of a minus sign (-).
This is a temporary syntax.

For example, the following command will list the available option for postprocessing:

.. code-block:: shell

   docker run -it -v $PWD:/localfs compass --postprocess /h

And this other one will postprocess result in ``output-vertical_column`` directory:

.. code-block:: shell

   docker run -it -v $PWD:/localfs compass --postprocess /s output-vertical_column

Advanced options
----------------

When mounting a volume you may experience problems if your host uid is not the user uid. In this case use the ``--compass-uid`` option.

.. code-block:: shell

   docker run -it -v /my/local/dir:/localfs compass --compass-uid $UID

Depending on your platform you may need to replace ``$UID`` with ``\`id -u $USER\```.

Adding the ``-p`` flag will run a parallel job with ``\`nproc\``` procs.

Running without internet / updating the local image
---------------------------------------------------

If you only want to dowload (or update) the image you can use the ``pull`` docker command (onced logged in to the gitlab registry):

.. code-block:: shell

   docker pull registry.gitlab.inria.fr/charms/compass:develop

`Docker_ComPASS_for_beginners.docx </uploads/235e455c7afaae43687ce36b033026fd/Docker_ComPASS_for_beginners.docx>`_

Parallel runs
-------------

Parallel runs with docker need a specific configuration of the devices of the host that the container can access.
Cf. the thread `here <https://github.com/open-mpi/ompi/issues/4948>`_.

A simple (but dangerous ?) workaround can be to give the container access to all devices on the host with the `--privileged <https://docs.docker.com/engine/reference/run/#runtime-privilege-and-linux-capabilities>`_ option.

Encapsulate docker call in a bash command
-----------------------------------------

To reduce typing work you can encapsulate calls to ``docker run`` inside a command. For example, on linux, you might have the following script somewhere on your file:

.. code-block:: bash

   #!/bin/bash
   current_directory=`pwd`
   docker run -it --privileged -v ${current_directory}:/localfs registry.gitlab.inria.fr/charms/compass:develop $@

Mapping volumes on Windows 7
----------------------------

On windows 7 don't forget that there is a 3 level russian doll: you run docker from your favorite console (cmd/power shell/git bash/...), to run a docker container on a virtual machine. So the mapping between file systems goes like that:

windows file system <-> virtual machine file system <-> container file system

If you used the virtualbox driver to create the virtual machine your ``C:\Users\username`` folder is usually mapped to the virtual machine folder ``/c/Users/username`` (beware of the case).

That is to say that if your user name is Toto and have a compass script whose path on Windows is ``C:\Users\Toto\path\to\my_wonderfull_script.py`` and assuming that ``compass`` refers to a valid container (obtnained by tagging an existing one) you can run compass on it with:

.. code-block:: shell

   docker run -it --volume /c/Users/Toto/path/to:/localfs compass my_wonderfull_script.py

.. _work:

Work environment
================

The *charms/compass/work-environment* container which is available in the project registry is designed to **compile and modify ComPASS and run simulations** through a container that ships all what is needed (compilation environment with dependencies and a few goodies such as scipy or matplotlib modules). Consequently it is a relatively big container but you will only need to install it once! The idea is to work on your local/host file system through the ``/localfs`` container volume (\ ``localfs``\ for *local file system*\ )

Let suppose that you work in a local directory called ``/my/local/dir``. If not already done clone the ComPASS git repository:

.. code-block:: shell

   cd /my/local/dir
   git clone git@gitlab.inria.fr:charms/ComPASS.git

or via https if you prefer:

.. code-block:: shell

   cd /my/local/dir
   git clone https://gitlab.inria.fr/charms/ComPASS.git

Then once you have installed docker (cf `previous section <#installation>`_\ ) you can pull the container from the registry (log in with you gitlab id and password):

.. code-block:: shell

   docker login registry.gitlab.inria.fr
   docker pull registry.gitlab.inria.fr/charms/compass/work-environment

Then you can possibly tag it to use shorter name:

.. code-block:: shell

   docker tag registry.gitlab.inria.fr/charms/compass/work-environment compass-shell

Container are created with a default user which is called *compass* and has uid 1000 (cf. miscellaneous/docker/base). When mounting a a volume you may experience problems if your host uid is not the user uid. In this case use the ``--user`` flag from the docker run command (cf. `docker documentation <https://docs.docker.com/engine/reference/run/#user>`_\ ). Now you can run the container interacting with the code that is on your local filesystem mouting ``/my/local/dir`` through the ``/localfs`` volume:

.. code-block:: shell

   docker run -it -v /my/local/dir:/localfs compass-shell --compass-uid $UID

Depending on your platform the variable ``UID`` may not be defined and you can replace it with ```id -u $USER` ``.

Check that you see everything that you have in ``/my/local/dir`` through the ``/localfs`` of the container, i.e. once in the container ``ls /localfs`` should list the content of ``/my/local/dir`` and esspecially the ``ComPASS`` source directory.

Now, on the container you can just compile and run ComPASS (let's suppose that you have 4 cores on you host you can use them on the container):

.. code-block:: shell

   cd /localfs
   mkdir build
   cd build
   cmake ../ComPASS
   make -j4 install
   cd ../ComPASS/test/bulk
   mpirun -n 4 python3 doublet_on_cartesian_grid.py

You can also use ``ccmake`` instead of ``cmake`` if you want to fine tune the compilation (select the physics and so on).

All changes that you make in the source files will remain on your local file system (and you should push them onto dedicated branches on the git repository! :smiley:). You will need to redo ``make install`` in the build directory, each time you reenter the container.
