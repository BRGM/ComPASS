Physical units
--------------

The
`International System of Units <https://en.wikipedia.org/wiki/International_System_of_Units>`_
is used throughout the code without specifying units.

For example:

   - distances are expressed in meters
   - permeabilities are expressed in :math:`m^2`
   - porosity has no units
   - times are expressed in seconds
   - ...

.. ifconfig:: versionlevel <= '4'

    To help users, some precomputed quantities are available in
    the `ComPASS.utils.units <https://github.com/BRGM/ComPASS/blob/v4.4.x/ComPASS/utils/units.py>`_
    module.

    For example, considering any function :code:`f(t)` that is expecting an argument in seconds,
    just import another duration from :code:`ComPASS.utils.units` to call :code:`f`:

    .. code-block:: python

        from ComPASS.utils.units import year
        f(10*year)

.. ifconfig:: versionlevel > '4'

    To help users, some precomputed quantities are available in
    the `compass_utils.units <https://gitlab.com/compass/compass-v5/compass-python-utils/-/blob/main/src/compass_utils/units.py?ref_type=heads>`_
    module.

    For example, considering any function :code:`f(t)` that is expecting an argument in seconds,
    just import another duration from :code:`compass_utils.units` to call :code:`f`:

    .. code-block:: python

        from compass_utils.units import year
        f(10*year)
