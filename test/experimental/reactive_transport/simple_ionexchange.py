#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

#
# The example in this file is from the paper
# AJ Valocchi, RL Street, PV Roberts Transport of ion-exchanging solutes in groundwater:
# Chromatographic theory and field simulation Water Resources Research 17 (5),
# 1517-1527, 1981
# See Fig. 1
#

import ComPASS
from ComPASS.utils.units import *
from ComPASS.timeloops import standard_loop
from ComPASS.timestep_management import FixedTimeStep, TimeStepManager
import ComPASS.timestep as timestep
from ComPASS.simulation_context import SimulationContext
from scipy.optimize import newton_krylov
import numpy as np

import sys
from copy import copy


Lx = 30
nx = 100

grid = ComPASS.Grid(
    shape = (nx, 1, 1),
    extent = (Lx, Lx/nx, Lx/nx),
)

muf = 1 # dynamic viscosity of pure water at 100°C
pleft, pright = 0.369 * Lx,  0  # Darcy ve;ocity given 0.369 m / h
k_reservoir = 1 # m^2
omega_reservoir = 0.2
vdarcy = (k_reservoir /muf) *(pleft - pright) / Lx
print('Darcy velocity: ', vdarcy , 'm / h')
volcell =   (Lx/nx)**3

rhow = 1
b = 1
rhocp = 1 # g / l
K_reservoir = 0.069  # m^2 / h
rhos = 2 # g /l ?
Keq = 1

# All in mmol / l, last in mmol / g
ctot_init = 100
ctot_inj = 10
c1_init = 20
c1_inj = 9.5
cTbar = 6

on_the_left = lambda x: x <= grid.origin[0]
on_the_right = lambda x: x >= grid.origin[0] + grid.extent[0]

def select_dirichlet_nodes():
    x = ComPASS.global_vertices()[:,0]
    return on_the_left(x) | on_the_right(x)

# concentrations (=temperature !) treated separately, as there are 2 species
def set_boundary_conditions():
    def set_states(states, x):
        left = on_the_left(x)
        states.p[left] = pleft
        right = on_the_right(x)
        states.p[right] = pright
        both = (left | right)
        states.context[both] = 1
        states.S[both] = 1
        states.C[both] = 1.

    set_states(ComPASS.dirichlet_node_states(), ComPASS.vertices()[:,0])

def set_states_inj(states, x, name):
    left = on_the_left(x)
    right = on_the_right(x)
    if name == 'cT':
        states.T[left] = ctot_inj
        states.T[right] = ctot_init
    elif name == 'c1':
        states.T[left] = c1_inj
        states.T[right] = c1_init
    else:
        raise SystemExit('wrong name %s % (name)') 
        
# concentrations (=temperature !) treated separately, as there are 2 species
def set_initial_values():
    def set_states(states, x):
        states.context[:] = 1
        states.p[:] = pleft + (pright - pleft) * (x - grid.origin[0]) / Lx
        states.S[:] = 1
        states.C[:] = 1.

    set_states(ComPASS.node_states(),  ComPASS.vertices()[:,0])
    set_states(ComPASS.cell_states(), ComPASS.compute_cell_centers()[:,0])


ComPASS.load_eos('linear_water')
ComPASS.set_output_directory_and_logfile(__file__)

fluid_properties = ComPASS.get_fluid_properties()
fluid_properties.specific_mass = rhow
fluid_properties.volumetric_heat_capacity = b
fluid_properties.dynamic_viscosity = muf
ComPASS.set_rock_volumetric_heat_capacity(rhocp)

ComPASS.init(
    mesh = grid,
    set_dirichlet_nodes = select_dirichlet_nodes,
    cell_porosity = omega_reservoir,
    cell_permeability = k_reservoir,
    cell_thermal_conductivity = K_reservoir,
)

set_initial_values()
set_boundary_conditions()

def retrieve_concentrations(): 
    # should work, no copy needed because of hstack
    return np.hstack((ComPASS.cell_states().T, ComPASS.node_states().T))

def set_concentrations(C):
    ComPASS.cell_states().T[:] = C[:nbCells]
    ComPASS.node_states().T[:] = C[nbCells:]

def set_source_term(S):
    cellVAGvolume = ComPASS.porovolfouriercell()
    nodeVAGvolume = ComPASS.porovolfouriernode()
    
    cellheatsource = ComPASS.cellthermalsource()
    cellheatsource[:] =  - cellVAGvolume /omega_reservoir * S[:nbCells]
    nodeheatsource = ComPASS.nodethermalsource()
    nodeheatsource[:] =  - nodeVAGvolume /omega_reservoir * S[nbCells:]

def clear_source_term():
    ComPASS.cellthermalsource()[:] = 0
    ComPASS.nodethermalsource()[:] = 0
    

def make_one_timestep(t, dt, cTold, c1old):

    ts_manager = FixedTimeStep(dt)
    #ts_manager = TimeStepManager(initial_timestep=720,)
    newton =  ComPASS.newton.Newton(1e-5, 20, ComPASS.newton.LinearSolver(1e-6, 150))  # ComPASS.default_Newton() #
    context = SimulationContext()
	
    # do cT first (linear)
    set_concentrations(cTold)
    set_states_inj(ComPASS.dirichlet_node_states(), ComPASS.vertices()[:,0], 'cT')
    clear_source_term()

    timestep.make_one_timestep(
        newton, ts_manager.steps(),
        simulation_context=context,
    )

    cTnew = retrieve_concentrations()

    print('---------- cT --> c1 ------------------')
    # Next do c1, solve non-linear system with Newton Krylov
    def freac(c1, cT):
        return Keq * rhos *cTbar *(c1 / (cT +(Keq-1)*c1))

    fc1old = freac(c1old, cTold)
    #fc1oldc = fc1old[:nbCells]
    #fc1oldn = fc1old[nbCells:]

    def fchim(cprev):
        set_concentrations(c1old)
        set_states_inj(ComPASS.dirichlet_node_states(), ComPASS.vertices()[:,0], 'c1')
        set_source_term((freac(cprev, cTnew) - fc1old) / dt)
        
        #np.set_printoptions(precision=8,linewidth=150)
        #print(heatsource)
        
        timestep.make_one_timestep(
            newton, ts_manager.steps(),
            simulation_context=context,
        )
        return retrieve_concentrations()
        
    cinit = c1old  #  np.random.rand(nx) 
    c1new = newton_krylov(lambda c: c-fchim(c), cinit, method='lgmres', verbose=1)

    #np.set_printoptions(precision=8,linewidth=150)
    #print(cTnew[:nbCells])
    #print(cTnew[nbCells:])
    
    return cTnew, c1new

##def plot_1D_concentrations(t, conc):
def plot_1D_concentrations(t, cT, c1):
    xc = ComPASS.compute_cell_centers()[:,0]
    xn = ComPASS.vertices()[:,0]
    ##conc = ComPASS.cell_states().T
    import ComPASS.utils.mpl_backends as mpl_backends
    plt = mpl_backends.import_pyplot(False)
    if plt:
        plt.clf()
        plt.subplot(221)
        plt.plot(xc, cT[:nbCells])
        plt.subplot(222)
        plt.plot(xn[:nbCells], cT[nbCells:2*nbCells])
        plt.xticks(np.arange(0,Lx+.01,step=5))
        plt.title('t='+str(t))
        plt.ylabel('Total conc')
        plt.subplot(223)
        plt.plot(xc, c1[:nbCells])
        plt.xticks(np.arange(0,Lx+.01,step=5))
        plt.subplot(224)
        plt.plot(xn[:nbCells], c1[nbCells:2*nbCells])
        plt.xticks(np.arange(0,Lx+.01,step=5))
        plt.xlabel('x in meters')
        plt.ylabel('C1 conc')
        plt.draw()
        plt.pause(0.1)

        plt.savefig(ComPASS.to_output_directory('conc_'+str(t)), format='png')

t=0
final_time = 55 # hrs
dt = 0.2

nbNodes = ComPASS.global_number_of_nodes()
nbCells = ComPASS.global_number_of_cells()
nbDofs = nbNodes + nbCells

cTold = np.tile(ctot_init, nbDofs)
c1old = np.tile(c1_init, nbDofs)


#midcurve= []
while t<final_time :
    print("===== Doing time ", t)

    cTnew, c1new = make_one_timestep (t, dt, cTold, c1old)
    plot_1D_concentrations(t, cTnew, c1new )
    t =t + dt
    cTold = cTnew
    c1old = c1new
#    ind = np.where(cTold['cells'] > 45)[0][0]
#    midcurve.append([ind, t, cTold['cells'][ind]])

#print(*midcurve, sep='\n')
