Documentation
=============

The documentation is generated using `Sphinx <https://www.sphinx-doc.org/>`_
and the lightweight markup language *reStructuredText* (cf. the `primer <https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html>`_).
Any document referenced directly or indirectly through the `docs/index.rst` file
will be included in the generated documentation and served on gitlab pages for the main development branch.

The documentation also includes documentation generated from python docstrings
using the sphinx `autodoc <https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html>`_ extension
and C++ and Fortran code using `doxygen <https://www.doxygen.nl/index.html>`_.

For test purpos, documentation can be generated using the `registry.gitlab.inria.fr/charms/compass/doc-environment`
docker image and executing the generate `generate_doc.bash` script from within the `docs` directory.
The generated documentation will output to the `docs/html` directory.

To generate the documentation locally, execute:

.. code:: bash

    cd doc
    cp -vf ../README.rst ../LICENSE.rst .
    cp -vf ../CodingConventions.rst ./developers_corner/
    sphinx-apidoc ../ComPASS -o python_reference
    sphinx-build . html
