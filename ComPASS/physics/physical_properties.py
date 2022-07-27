#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import numpy as np
from numba import cfunc, jit, carray, vectorize
from numba.types import void, double, CPointer as p, Record, NestedArray
from typing import NamedTuple
from types import FunctionType


class PhaseProperty(NamedTuple):
    with_derivatives: FunctionType
    without_derivatives: FunctionType


# one by phase, n_comp can differ ?
class PhaseStateStruct:
    def __init__(self, n_comp):
        self.n_comp = n_comp
        self.Xt = Record.make_c_struct(
            [
                ("pressure", double),
                ("temperature", double),
                (
                    "molar_fractions",
                    NestedArray(dtype=double, shape=(n_comp,)),
                ),
            ]
        )

    def Xalpha(self, pressure, temperature, molar_fractions=None):
        if molar_fractions is None:
            assert self.n_comp == 1, "Molar fractions missing in Xalpha"
            molar_fractions = np.ones(self.n_comp)
        if (
            np.isscalar(pressure)
            and np.isscalar(temperature)
            and np.ndim(molar_fractions) == 1
        ):
            return np.rec.array((pressure, temperature, molar_fractions), dtype=self.Xt)
        else:  # vector
            array_size = max(np.size(pressure), np.size(temperature))
            if np.ndim(molar_fractions) > 1:
                array_size = max(array_size, np.shape(molar_fractions)[0])
            X = self.empty_Xalpha(array_size)
            X["pressure"] = pressure
            X["temperature"] = temperature
            X["molar_fractions"] = molar_fractions
            return X

    def empty_Xalpha(self, len=None):
        if len is not None and len > 1:
            return np.empty(len, dtype=self.Xt).view(np.recarray)
        return np.empty((), dtype=self.Xt).view(np.recarray)


# contains lists with the phases properties
class FluidMixtureProperties:

    __property_names__ = ("dynamic_viscosity", "molar_density")

    def __init__(self, n_phases, n_components):
        self.n_components = n_components
        self.phase_state_type = PhaseStateStruct(n_components)
        self.properties = {
            name: [None] * n_phases  # list of CompiledPhaseProperty, one by phase
            for name in self.__property_names__
        }

    def register(self, phase_property_name, phase_id, phase_property):
        self.properties[phase_property_name][phase_id] = CompiledPhaseProperty(
            self.phase_state_type,
            phase_property.with_derivatives,
            phase_property.without_derivatives,
        )

    def __getattr__(self, name):
        if name in self.__property_names__:
            return tuple(self.properties[name])
        raise AttributeError(name)


class CompiledPhaseProperty:
    def __init__(
        self,
        phase_state_type,
        phase_law_with_derivatives,
        phase_law_without_derivatives,
    ):
        # the functions are stored so the address always exist
        self.phase_state_type = phase_state_type
        self.py_with_derivatives = phase_law_with_derivatives
        self.py_without_derivatives = phase_law_without_derivatives
        # The compilation is done when calling for the first time
        # the compiled function because the python function can be
        # set several times (init with the default one then the user
        # can give an other one)
        self._c_with_derivatives_compiled = None
        self._c_without_derivatives_compiled = None
        self._vect_with_derivatives_compiled = None
        self._vect_without_derivatives_compiled = None

    @property
    def c_with_derivatives(self):
        if self._c_with_derivatives_compiled is None:
            self._c_with_derivatives_compiled = self._compile_c_with_derivatives()
        return self._c_with_derivatives_compiled

    @property
    def c_without_derivatives(self):
        if self._c_without_derivatives_compiled is None:
            self._c_without_derivatives_compiled = self._compile_c_without_derivatives()
        return self._c_without_derivatives_compiled

    @property
    def vect_with_derivatives(self):
        if self._vect_with_derivatives_compiled is None:
            self._vect_with_derivatives_compiled = self._compile_vect_with_derivatives()
        return self._vect_with_derivatives_compiled

    @property
    def vect_without_derivatives(self):
        if self._vect_without_derivatives_compiled is None:
            self._vect_without_derivatives_compiled = (
                self._compile_vect_without_derivatives()
            )
        return self._vect_without_derivatives_compiled

    def _compile_c_without_derivatives(self):
        function_signature = double(p(self.phase_state_type.Xt))  # val = func(Xalpha)
        make_func_cfunc = cfunc(function_signature, nopython=True, cache=True)
        c_f = jit(
            self.py_without_derivatives, nopython=True
        )  # nopython=True by default (soon)

        def c_wrapped_function(X_):
            X = carray(X_, 1)[0]
            return c_f(X)

        return make_func_cfunc(c_wrapped_function)

    def _compile_c_with_derivatives(self):
        subroutine_signature = void(
            p(self.phase_state_type.Xt),  # Xalpha
            p(double),  # f
            p(self.phase_state_type.Xt),  # dfdXalpha
        )
        make_sub_cfunc = cfunc(subroutine_signature, nopython=True, cache=True)
        c_f = jit(self.py_with_derivatives, nopython=True)

        def c_wrapped_subroutine(X_, val_, dfdX_):
            X = carray(X_, 1)[0]
            val = carray(val_, 1)
            dfdX = carray(dfdX_, 1)[0]
            val[0] = c_f(X, dfdX)

        return make_sub_cfunc(c_wrapped_subroutine)

    def _compile_vect_with_derivatives(self):
        return vectorize([double(self.phase_state_type.Xt, self.phase_state_type.Xt)])(
            self.py_with_derivatives
        )

    def _compile_vect_without_derivatives(self):
        return vectorize([double(self.phase_state_type.Xt)])(
            self.py_without_derivatives
        )
