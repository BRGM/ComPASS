#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import ComPASS
import importlib
#import doublet_utils
from ComPASS.utils.units import *
import ComPASS.timestep as timestep
from  ComPASS.timeloops import standard_loop
from ComPASS.timestep_management import FixedTimeStep, TimeStepManager
from ComPASS.simulation_context import SimulationContext
from ComPASS.newton import Newton, LinearSolver
import ComPASS.mpi as mpi
import matplotlib.pyplot as plt
import scipy.sparse as sps
import numpy as np


class Transport (object) :
	def __init__(self, rhow, b, rhocp, muf, porosity, perm, conduct, onecomp, exact_sol, grid, pleft, pright, a, CL, CR):
	
		if onecomp:
		    if exact_sol:
		    	ComPASS.load_eos('linear_water')
		    else:
		    	ComPASS.load_eos('liquid_water')
		else:
		    ComPASS.load_eos('water_with_tracer')
		    
		fluid_properties = ComPASS.get_fluid_properties()
		fluid_properties.specific_mass = rhow
		fluid_properties.volumetric_heat_capacity = b
		fluid_properties.dynamic_viscosity = muf
		ComPASS.set_rock_volumetric_heat_capacity(rhocp)

		self.on_the_left = lambda x: x <= grid.origin[0]
		self.on_the_left_in_the_middle = lambda x,y: (x <= grid.origin[0] ) & (grid.extent[1]/2-a <y) & (y< grid.extent[1]/2+a)
		self.on_the_right = lambda x: x >= grid.origin[0] + grid.extent[0]

		# only sequential allowed here
		assert ComPASS.mpi.master_proc_rank == ComPASS.mpi.proc_rank
		ComPASS.set_output_directory_and_logfile(__file__)
		
		ComPASS.init(
		    grid = grid,
		    set_dirichlet_nodes = self.select_dirichlet_nodes,
		    cell_porosity = porosity,
		    cell_permeability = perm,
		    cell_thermal_conductivity = conduct,
		)

		self.nb_nodes = ComPASS.global_number_of_nodes()
		self.nb_cells = ComPASS.global_number_of_cells()
		self.nb_points = self.nb_nodes + self.nb_cells
		self.CL = CL
		self.CR = CR
		self.grid =grid
				
		self.set_initial_values(pright, onecomp)
		self.set_boundary_conditions(pleft, pright, onecomp)

	def select_dirichlet_nodes(self):
	    x = ComPASS.global_vertices()[:,0]
	    return self.on_the_left(x) | self.on_the_right(x)
	
	def set_boundary_conditions(self, pleft, pright, onecomp):
	    def set_states(states, x,y):
	    	left = self.on_the_left(x)
	    	right = self.on_the_right(x)
	    	leftmid = self.on_the_left_in_the_middle(x, y)
	    	both = (left | right)

	    	states.p[left] = pleft
	    	states.p[right] = pright

	    	states.context[both] = 1
	    	states.S[both] = 1
	    	if onecomp:
	    		states.C[both] = 1.
	    	else:
	    		states.C[left] = (1, 0) #(0., 1)
	    		states.C[leftmid] = (0., 1.)
	    		states.C[right] = (1, 0)
	    verts = ComPASS.vertices()
	    set_states(ComPASS.dirichlet_node_states(), verts[:,0], verts[:,1])

	def set_states_inj(self, states, x,y, idx):
		left = self.on_the_left(x)
		right = self.on_the_right(x)
		leftmid = self.on_the_left_in_the_middle(x, y)
		both = (left | right)

		states.T[right] = self.CR[idx]
		states.T[left] =  self.CR[idx]
		states.T[leftmid]  = self.CL[idx] 

	def set_initial_values(self, p0, onecomp):
	    def set_states(states, x):
	    	states.context[:] = 1
	    	states.p[:] = p0 # pleft + (pright - pleft) * (x - grid.origin[0]) / Lx
	    	states.T[:] = 0 #Tright
	    	states.S[:] = 1
	    	if onecomp:
	    		states.C[:] = 1.
	    	else:
	    		states.C[:] = (1, 0)
	     
	    verts = ComPASS.vertices()
	    cellcenters = ComPASS.compute_cell_centers()
	    set_states(ComPASS.node_states(),  verts[:,0])
	    set_states(ComPASS.cell_states(), cellcenters[:,0])

	def plot_concentrations(self, t, conc) :

	    Caq=["Ca2+", "Na+", "Cl-", "K+"]
	    xy = ComPASS.compute_cell_centers()[:,0:2]
	    nx = self.grid.shape[0] 
	    ny = self.grid.shape[1] 
	    XX = xy[:, 0].reshape(ny, nx)
	    YY = xy[:, 1].reshape(ny, nx)
	    fig = plt.figure(1)
	    
	    for i in range (len(Caq)) :
		    Ci_nodes = conc[i][0:self.nb_nodes] 
		    Ci_cells = conc[i][self.nb_nodes : self.nb_points] 

		    ax = plt.subplot(4, 1, i+1)
		    #lvs = [1e-10, 0.000147, 0.000294, 0.000442, 0.000589]
		    #cs =plt.contourf(XX,YY,np.reshape(Ci_cells,[ny,nx]), 15, levels=lvs, vmin=1e-10, vmax=0.000589,cmap='jet')
		    cs =plt.contourf(XX,YY,np.reshape(Ci_cells,[ny,nx]), 15, cmap='jet')
		    fig.colorbar(cs)  
		    plt.setp(ax.get_xticklabels(), visible=False)
		    plt.gcf().subplots_adjust(hspace = 0.4)
		    plt.title(Caq[i])
		    plt.suptitle('t='+str(t/3600.)+ '  (hrs)')
	    plt.setp(ax.get_xticklabels(), visible=True)
	    plt.xlabel('x in meters')
	    plt.ylabel('y in meters') 
	    plt.draw()
	    plt.pause(0.1)
	    #if ((t==0)|(t==21600)|(t==28800)) :
	    #	plt.savefig(ComPASS.to_output_directory('diff_cell_conc_at_time_'+str(t)),format='png')
	    plt.savefig(ComPASS.to_output_directory('cell_conc_at_time_'+str(t)),format='png')
	    plt.clf()
	    
	def plot_1D_concentrations(self, t, conc) :

	    Caq=["Ca2+", "Na+", "Cl-", "K+"]
	    x = ComPASS.cell_centers()[:, 0]
	    nx = self.grid.shape[0] 
	    ny = self.grid.shape[1] 
	    ligne = nx*int(ny/2)
	    fig2 = plt.figure(2)
	    
	    for i in range (len(Caq)) :
		    Ci_nodes = conc[i][0:self.nb_nodes] 
		    Ci_cells = conc[i][self.nb_nodes : self.nb_points] 

		    ax = plt.subplot(2, 1, 1)
		    #cs =plt.plot(x[ligne:ligne+30], Ci_cells[ligne:ligne+30], label=Caq[i])
		    cs =plt.plot(x[ligne:ligne+nx], Ci_nodes[ligne:ligne+nx], label=Caq[i])
		    plt.title('plot at nodes')
		    plt.suptitle('t='+str(t/3600.)+ ' (hrs)')

		    ax2 = plt.subplot(2, 1, 2)
		    cs =plt.plot(x[ligne:ligne+nx], Ci_cells[ligne:ligne+nx], label=Caq[i])
		    plt.title('plot at cells')
		    plt.title('t='+str(t/3600.)+ ' (hrs)')
	    plt.xlabel('x in meters')
	    plt.ylabel('conc') 
	    plt.legend()
	    plt.draw()
	    plt.pause(0.1)
	    plt.savefig(ComPASS.to_output_directory('1D_cell_conc_at_time_'+str(t)),format='png')
	    #if ((t==0)|(t==21600)|(t==28800)) :
	    #	plt.savefig(ComPASS.to_output_directory('diff_1D_cell_conc_at_time_'+str(t)),format='png')
	    plt.clf()

	# the function that advects and diffuses all concentration
	def transport_concentrations(self, t, ts_manager, Cold, srcF):
	    Nc = Cold.shape[0]  
	    Cnew = np.zeros_like(Cold)
	    
	    compute_all_concentrations = True

	    newton =  ComPASS.default_Newton() #Newton(1e-5, 5, LinearSolver(1e-6, 150))  #
	    context = SimulationContext()
		
	    # while loop to find suitable dt (other strategies are possible...)
	    while compute_all_concentrations:

		    for i in range(Nc) :
		    	print ("espece ========> ", i)
		    	#print ("ts_manager.current_step ========> ", ts_manager.current_step)

		    	Src_nodes =  srcF[i][0:self.nb_nodes]
		    	Src_cells =  srcF[i][self.nb_nodes : self.nb_points] 
		    
		    	if mpi.is_on_master_proc:
	    		
		    		#buffer1 = np.array(getattr(ComPASS.kernel, 'get_node_heat_source_buffer')(), copy = False)
		    		buffer1 = ComPASS.nodethermalsource()
		    		buffer1[:] = Src_nodes/ts_manager.current_step 

		    		#buffer2 = np.array(getattr(ComPASS.kernel, 'get_cell_heat_source_buffer')(), copy = False)
		    		buffer2 = ComPASS.cellthermalsource()
		    		buffer2[:] = Src_cells/ts_manager.current_step

		    	ComPASS.node_states().T[:] = Cold[i][0:self.nb_nodes]
		    	ComPASS.cell_states().T[:] = Cold[i][self.nb_nodes : self.nb_points]   

		    	self.set_states_inj(ComPASS.dirichlet_node_states(), ComPASS.vertices()[:,0], ComPASS.vertices()[:,1], i)

		    	deltat = timestep.make_one_timestep(
		            newton, ts_manager.steps(),
		            simulation_context=context,
		    	)
		                        	       	
		    	Cnew[i, 0:self.nb_nodes] = ComPASS.node_states().T
		    	Cnew[i, self.nb_nodes : self.nb_points] = ComPASS.cell_states().T		    	
		    	
		    else:
		    	compute_all_concentrations = False #job is done all concentrations have been computed with dt
	    #self.plot_1D_concentrations(t, Cnew)
		    	
	    return Cnew

