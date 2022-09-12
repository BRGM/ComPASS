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
  \text{set of phases}      &   \mathcal{P}   &   \emptyset    &   \\
  \hdashline
  \text{set of components}  &   \mathcal{C}   &   \emptyset    &   \\
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
  \text{permeability tensor} & \Lambda & m^2 &  \\
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
  \text{partial molar enthalpy} & h_i^\alpha, i \in \mathcal{C}, \alpha \in \mathcal{P}   & J.mol^{-1} & h_i^\alpha = \left.\frac{\partial H^\alpha}{\partial N_i}\right|_{T,P}  \text{with } H^\alpha \\
   & & & \text{ the extensif enthalpy in } J \\
  \hdashline
  \text{molar enthalpy} & h^\alpha, \alpha \in \mathcal{P}   & J.mol^{-1} & h^\alpha = \sum\limits_{i\in\mathcal{C}} C_i^\alpha h_i^\alpha \\
  \hdashline
  \text{molar internal energy} & e^\alpha, \alpha \in \mathcal{P}   & J.mol^{-1} &  \\
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
