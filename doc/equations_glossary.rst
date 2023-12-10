Physical equations solved by ComPASS
====================================

Balance equations
-----------------

ComPASS solves the following molar and thermal balance equations:

.. math::

  \begin{array}{r l c}
  \phi \partial_t n_i + \mathrm{div} ( \sum\limits_{\alpha\in \mathcal{P}} C^\alpha_i \zeta^\alpha {\mathbf q}^\alpha ) &=& 0, & i\in \mathcal{C}, \\
  \phi \partial_t E_f + (1-\phi) \partial_t E_r + \mathrm{div} ( \sum\limits_{\alpha \in \mathcal{P}} h^\alpha \zeta^\alpha {\mathbf q}^\alpha - \lambda \nabla T ) &=& 0, &
  \end{array}

where

.. math::

  {\mathbf q}^\alpha = -{kr^\alpha \over \mu^\alpha } \Lambda (\nabla P^\alpha - \rho^\alpha {\mathbf g}).

The previous system of conservation equations is coupled to the following local closure laws:

.. math::

  \begin{array}{r l c}
  & \sum\limits_{\alpha\in \mathcal{P}} S^\alpha = 1, \\
  & \sum\limits_{i \in \mathcal{C}} C^{\alpha}_i = 1, \quad \alpha\in\mathcal{P}, \\
  & P^\alpha - P^\beta = P_c^{\alpha,\beta}, \quad \alpha\neq \beta, \quad (\alpha,\beta) \in \mathcal{P}^2, \\
  & f_i^\alpha = f_i^\beta, \quad \alpha \neq \beta, \quad (\alpha, \beta) \in \mathcal{P}^2.
  \end{array}

Glossary and notations
----------------------

.. math::

  \begin{array}{| c | c | c |}
  \hline
  \text{name}                  &     \text{notation}       &   \text{unit} & \text{comments} \\
  \hline
  \text{set of phases}      &   \mathcal{P}   &   \emptyset    & \text{liquid and gas for example} \\
  \hdashline
  \text{set of components}  &   \mathcal{C}   &   \emptyset    & \text{water and air for example} \\
  \hdashline
  \text{porosity}           &   \phi          &   \emptyset    & \text{void fraction} \\
   & & & \text{of the porous medium}  \\
  \hdashline
  \text{molar mass}    & M_i, i \in \mathcal{C} &   kg.mol^{-1}    &  \text{mass of one mole of } i   \\
  \hdashline
  \text{saturation}           &   S^\alpha &   \emptyset    & \text{phase fraction:} \sum\limits_{\alpha \in \mathcal{P}} S^\alpha = 1  \\
  \hdashline
  \text{molar fractions} & C_i^\alpha, i \in \mathcal{C}, \alpha \in \mathcal{P} & \emptyset & \sum\limits_{i\in\mathcal{C}} C_i^\alpha = 1, \alpha \in \mathcal{P} \\
  \hdashline
  \text{number of moles}    & n_i, i \in \mathcal{C} &   mol.m^{-3}    &  \text{number of moles}  \\
   &  &       &   \text{per unit pore volume}  \\
   &  &       &   n_i =  \sum\limits_{\alpha \in \mathcal{P}} S^\alpha C_i^\alpha \zeta^\alpha \\
  \hdashline
  \text{temperature} & T & K &  \\
  \hdashline
  \text{pressure} & P^\alpha, \alpha \in \mathcal{P} & Pa &  \\
  \hdashline
  \text{gravity} & {\mathbf g}   & m.s^{-2} &  \\
  \hdashline
  \text{bulk thermal} & \lambda & W.m^{-1}.K^{-1} & \text{of the fluid} \\
  \text{conductivity} & & & \text{and rock mixture} \\
  \hdashline
  \text{permeability tensor} & \Lambda & m^2 & \text{resistance of the} \\
    & & & \text{medium to the flow}  \\
  \hdashline
  \text{relative permeability} & kr^\alpha, \alpha \in \mathcal{P} & \emptyset &  \\
  \hdashline
  \text{dynamic viscosity} & \mu^\alpha, \alpha \in \mathcal{P} & Pa.s &  \\
  \hdashline
  \text{capillary pressure} & P_c^{\alpha,\beta}, (\alpha,\beta) \in \mathcal{P}^2 & Pa &  \\
  \hdashline
  \text{molar density} & \zeta^\alpha, \alpha \in \mathcal{P}   & mol.m^{-3} & \text{number of moles} \\
   & & & \text{per unit pore volume} \\
  \hdashline
  \text{volumetric mass} & \rho^\alpha, \alpha \in \mathcal{P}   & kg.m^{-3} & \text{mass per unit pore volume} \\
  \text{density} &  &       &   \rho^\alpha = (\sum\limits_{i\in\mathcal{C}} C_i^\alpha  M_i) \zeta^\alpha  \\
  \hdashline
  \text{fugacity} & f_i^\alpha, i \in \mathcal{C}, \alpha \in \mathcal{P}  & Pa &  \\
  \hdashline
  \text{specific enthalpy} & H^\alpha, \alpha \in \mathcal{P}   & J.kg^{-1} &  \\
  \hdashline
  \text{partial molar enthalpy} & h_i^\alpha, i \in \mathcal{C}, \alpha \in \mathcal{P}   & J.mol^{-1} & h_i^\alpha = \left.\frac{\partial \tilde{H}^\alpha}{\partial N_i}\right|_{T,P}  \text{with } \tilde{H}^\alpha \\
   & & & \text{ the extensif enthalpy in } J \\
  \hdashline
  \text{molar enthalpy} & h^\alpha, \alpha \in \mathcal{P}   & J.mol^{-1} & h^\alpha = \sum\limits_{i\in\mathcal{C}} C_i^\alpha h_i^\alpha \\
   & & & \text{and } h^\alpha = (\sum\limits_{i\in\mathcal{C}} M_i C_i^\alpha) H^\alpha \\
  \hdashline
  \text{internal energy} & E^\alpha, \alpha \in \mathcal{P}   & J.kg^{-1} &  \\
  \hdashline
  \text{molar internal energy} & e^\alpha, \alpha \in \mathcal{P}   & J.mol^{-1} & e^\alpha = (\sum\limits_{i\in\mathcal{C}} M_i C_i^\alpha) E^\alpha \\
  \hdashline
  \text{rock energy} & E_r  & J.m^{-3} &  \\
  \text{per unit rock volume} & & & \\
  \hdashline
  \text{fluid energy} & E_f  & J.m^{-3} & E_f = \sum\limits_{\alpha \in \mathcal{P}} S^\alpha \zeta^\alpha e^\alpha \\
  \text{per unit pore volume} & & & \\
  \hdashline
  \text{generalized Darcy} & {\mathbf q}^\alpha, \alpha \in \mathcal{P}   & m.s^{-1} &  \\
  \text{velocity} & &  &  \\
  \hdashline
  \text{molar flux} & q_i, i \in \mathcal{C} & mol.m^{-2}.s^{-1} & q_i = \sum\limits_{\alpha \in \mathcal{P}} C_i^\alpha \zeta^\alpha {\mathbf q}^\alpha\\
  \hdashline
  \text{energy flux} & q_e   & J.m^{-2}.s^{-1} & q_e = \sum\limits_{\alpha \in \mathcal{P}} h^\alpha \zeta^\alpha {\mathbf q}^\alpha  \\
  \hline
  \end{array}

You can find information about the default physical properties in the :ref:`physics section<Available physics>`.
The default values depends on which physics you use.

.. _water2ph_equations:

.. warning::

  In some physics which contain a **single component** (like *water2ph*) the balance
  equations are adapted to solve the mass and energy balance equations (by multiplying the first equation
  by the molar mass of the component, :math:`M_{H_2O}` for *water2ph*):


  .. math::

    \begin{array}{r l c}
    \phi \partial_t (\sum\limits_{\alpha\in P} \rho^\alpha S^\alpha C_i^\alpha) + \mathrm{div} ( \sum\limits_{\alpha\in \mathcal{P}} C^\alpha_i \rho^\alpha {\mathbf q}^\alpha ) &=& 0, & i\in \mathcal{C}, \\
    \phi \partial_t (\sum\limits_{\alpha\in P} \rho^\alpha S^\alpha E^\alpha) + (1-\phi) \partial_t E_r + \mathrm{div} ( \sum\limits_{\alpha \in \mathcal{P}} H^\alpha \rho^\alpha {\mathbf q}^\alpha - \lambda \nabla T ) &=& 0. &
    \end{array}

  The energy balance equation remains the same because :math:`\rho^\alpha E^\alpha = \zeta^\alpha e^\alpha`
  and :math:`H^\alpha \rho^\alpha = h^\alpha \zeta^\alpha`.

  To do so, the core of ComPASS remains identical and we define :math:`\zeta^\alpha = \rho^\alpha`, :math:`e^\alpha = E^\alpha`
  and :math:`h^\alpha = H^\alpha`, which is equivalent to consider that :math:`M_{H_2O}=1`
  (let us recall that :math:`C_{H_2O}^\alpha = 1` because there is only one component).

  .. ifconfig:: versionlevel <= '4'

    **Careful**: this can have an impact on the set-up of the simulation,
    especially when setting a
    :ref:`Neumann boundary flux<Neumann faces>`.
    Keep in mind that, when using the *water2ph* physics, you need to give
    the mass flux instead of the molar flux (:math:`M_{H_2O}=1`)
    using the :code:`ComPASS.NeumannBC().molar_flux` object.

    .. code-block:: python

        Neumann = ComPASS.NeumannBC()
        Neumann.heat_flux = bottom_heat_flux # in W/m^2 = J/m^2/s
        Neumann.molar_flux[:] = Qm # one value by component in kg/m^2/s !
        face_centers = simulation.face_centers()
        simulation.set_Neumann_faces(face_centers[:, 2] <= -H, Neumann)

  .. ifconfig:: versionlevel > '4'

    **Careful**: this can have an impact on the set-up of
    the simulation, especially when setting a
    :ref:`Neumann boundary flux<Neumann flux>`.
    Keep in mind that, when using the *water2ph* physics, you need to give
    the mass flux instead of the molar flux (:math:`M_{H_2O}=1`).

    .. code-block:: python

      model = Coats("water2ph")
      data = model.new_data(geom.mesh.dimension)
      Neumann_flux = data.boundary_conditions.Neumann_flux.dtype()
      Neumann_flux.energy = bottom_heat_flux  # in W/m^2 = J/m^2/s
      Neumann_flux.molar[model.components["water"]] = qmass  # in kg/m^2/s !
      bottom_faces = geom_traits.bottom_boundary_faces(geom)
      # init all the bottom faces with the same Neumann flux
      data.boundary_conditions.Neumann_flux[bottom_faces] = Neumann_flux

Available physics
-----------------

.. ifconfig:: versionlevel <= '4'

  .. include:: physics_summary_v4.rst

.. ifconfig:: versionlevel > '4'

  .. include:: physics_summary_v5.rst

For instructions to load the physics,
refer to :ref:`this section<Load the physics>`.

After loading a physics, you can also
:ref:`set your own physical properties<Fluid physical properties>`.

.. ifconfig:: versionlevel <= '4'

  .. include:: linear_water.rst

.. ifconfig:: versionlevel > '4'

  .. include:: pure.rst

water2ph
........

This physics contains two phases (by default liquid and gas),
one component (by default water).

.. warning::
  In this physics, to solve the mass balance equation instead of the
  molar balance equation, the molar and volumetric mass densities are considered equal,
  and the water molar mass is set to 1.
  For more information, refer to :ref:`this paragraph<water2ph_equations>`.

The default physical properties are

* gas molar and volumetric mass densities:

.. math::

    \begin{array}{ r l }
        \zeta^g =& \rho^g  \\
        M_{H_2O} =& 0.018016  \\
        R =& 8.3145  \\
        Z =& 1  \\
        \rho^g =& {\mathbf{P^g} * M_{H_2O} \over{R * \mathbf{T} * Z}}
    \end{array}

* liquid molar and volumetric mass densities:

.. math::

    \begin{array}{ r l }
      \zeta^l =& \rho^l  \\
      rho_0 =& 780.83795  \\
      a =& 1.6269192  \\
      b =& -3.0635410e^{-3}  \\
      a1 =& 2.4638e^{-9}  \\
      a2 =& 1.1343e^{-17}  \\
      b1 =& -1.2171e^{-11}  \\
      b2 =& 4.8695e^{-20}  \\
      c1 =& 1.8452e^{-14}  \\
      c2 =& -5.9978e^{-23}  \\
      T_{square} =& \mathbf{T}^2  \\
      ss =& rho_0 + a * \mathbf{T} + b * T_{square}  \\
      cw =& (
          a1
          + a2 * \mathbf{P^l}
          + \mathbf{T} * (b1 + b2 * \mathbf{P^l})
          + T_{square} * (c1 + c2 * \mathbf{P^l})
      )  \\
      psat =& 1e^{-3} * (\mathbf{T} - 273)^4  \\
      p_{rel} =& \mathbf{P^l} - psat  \\
      \rho^l =& ss * (1 + cw * p_{rel})  \\
    \end{array}

* components molar masses:

.. math::

  M_{H_2O} = 1

* gas viscosity:

.. math::

  \mu^g = (0.361 * \mathbf{T} - 10.2) * 1.0e^{-7}

* liquid viscosity:

.. math::

    \begin{array}{ r l }
        Tref =& \mathbf{T} - 273 - 8.435  \\
        b =& \sqrt{8078.4 + Tref^2}  \\
        a =& 0.021482 * (Tref + b) - 1.2  \\
        \mu^l =& 1.0e^{-3} / a  \\
    \end{array}

* gas molar enthalpy:

.. math::

    h^g = 1990.89e^3 + 190.16e^1 * \mathbf{T}


* liquid molar enthalpy (pure liquid water molar enthalpy):

.. math::

    \begin{array}{ r l }
        a =& -14.4319e^3 \\
        b =& 4.70915e^3 \\
        cc =& -4.87534 \\
        d =& 1.45008e^{-2} \\
        T_0 =& 273 \\
        TdegC =& \mathbf{T} - T_0 \\
        ss =& a + b * TdegC + cc * TdegC^{2} + d * TdegC^{3} \\
        h^l =& M_{H_2O} * ss
    \end{array}

* saturation pressure:

.. math::

    psat = 1e^{-3} * (\mathbf{T} - 273)^4

* relative permeabilities:

.. math::

    kr^{\alpha} = \mathbf{S^\alpha}^2

* the capillary pressure is null:

.. math::

  \mathbf{P^g} = \mathbf{P^l}

* the rock volumetric heat capacity is :math:`1.6e^6` (:math:`800 * 2000 \,\, J/m^3`).


immiscible2ph
.............

This physics contains two phases (by default liquid and gas),
two components (by default water and air), only water in liquid phase
and air in gas phase.

.. ifconfig:: versionlevel > '4'

  This physics is not implemented yet in ComPASS v5.


.. The default physical properties are


diphasic
........

This physics contains two phases (by default liquid and gas), two components
(by default water and air), all components can be in all phases.

The default physical properties are

* molar densities:

.. math::

    \begin{array}{ r l }
        \zeta^l =& 1000 / 18e^{-3} \\
        \zeta^g =& {\mathbf{P^g} \over{Rgp * \mathbf{T}}}
    \end{array}

* components molar masses:

.. math::

    \begin{array}{ r l }
        M_{H_2O} =& 18e^{-3} \\
        M_{air} =& 29e^{-3}
    \end{array}

* viscosities:

.. math::

    \begin{array}{ r l }
        \mu^l =& 1e^{-3} \\
        \mu^g =& 2e^{-5}
    \end{array}

* gas molar enthalpy:

.. math::

    \begin{array}{ r l }
        a =& 1990.89e^3 \\
        b =& 190.16e^3 \\
        cc =& -1.91264e^3 \\
        d =& 0.2997e^3 \\
        Cp^g =& 1000 \\
        T_s =& \mathbf{T} / 100 \\
        ss =& a + b * T_s + cc * T_s^{2} + d * T_s^{3} \\
        \beta_{air} =& Cp^g * M_{air} * ss \\
        \beta_{H_2O} =& M_{H_2O} * \mathbf{T} \\
        h^g =& \sum\limits_{i\in\mathcal{C}} \beta_i C_i^\alpha
    \end{array}


* liquid molar enthalpy (pure liquid water molar enthalpy):

.. math::

    \begin{array}{ r l }
        a =& -14.4319e^3 \\
        b =& 4.70915e^3 \\
        cc =& -4.87534 \\
        d =& 1.45008e^{-2} \\
        T_0 =& 273 \\
        TdegC =& \mathbf{T} - T_0 \\
        ss =& a + b * TdegC + cc * TdegC^{2} + d * TdegC^{3} \\
        h^l =& M_{H_2O} * ss
    \end{array}

* saturation pressure:

.. math::

    psat = 100 * \exp{(46.784 - 6435 / \mathbf{T} - 3.868 * \log(\mathbf{T}))}

* relative permeabilities:

.. math::

    kr^{\alpha} = \mathbf{S^\alpha}^2

* the capillary pressure is null:

.. math::

  \mathbf{P^g} = \mathbf{P^l}

* the rock volumetric heat capacity is :math:`1.6e^6` (:math:`800 * 2000 \,\, J/m^3`).
