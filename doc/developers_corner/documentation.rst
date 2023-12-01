Documentation
=============

The documentation is generated using `Sphinx <https://www.sphinx-doc.org/>`_
and the lightweight markup language *reStructuredText*
(cf. the `primer <https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html>`_).
Any document referenced directly or indirectly through the `doc/index.rst` file
will be included in the generated documentation and served on gitlab pages
for the main development branch.

The documentation also includes documentation generated from python docstrings
using the sphinx `autodoc <https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html>`_
extension and C++ and Fortran code using `doxygen <https://www.doxygen.nl/index.html>`_.

For test purpose, documentation can be generated using the
`registry.gitlab.inria.fr/charms/compass/doc-environment`
docker image and executing the generate `make_the_doc` script
from within the `doc` directory:

.. ifconfig:: versionlevel <= '4'

    .. code:: bash

        source doc/make_the_doc --with-doxygen

.. ifconfig:: versionlevel > '4'

    .. code:: bash

        source doc/make_the_doc --with-doxygen --version 5

The generated documentation will output to the `doc/html` directory.
