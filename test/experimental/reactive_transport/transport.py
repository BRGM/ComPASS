#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import ComPASS
import importlib

# import doublet_utils
from ComPASS.utils.units import *
import ComPASS.timestep as timestep
from ComPASS.timeloops import standard_loop
from ComPASS.timestep_management import FixedTimeStep, TimeStepManager
from ComPASS.simulation_context import SimulationContext
from ComPASS.newton import Newton
from ComPASS.linalg import linear_solver
import ComPASS.mpi as mpi
import matplotlib.pyplot as plt
import scipy.sparse as sps
import numpy as np

# from mayavi import mlab


"""def variable_conductivity():
    xyz = self.simulation.compute_global_cell_centers()
    nbcells = xyz.shape[0]
    x = xyz[:, 0]
    Kleft = 5.55556e-9
    Kright = 5.55556e-9
    xmin = x.min()
    xmax = x.max()
    Kscalar = (x - xmin) * (Kright / (xmax - xmin)) + (xmax - x) * (
                                                                    Kleft / (xmax - xmin)
                                                                    )
    return np.array([np.diag([Ki, 0.1 * Ki, 0.1 * Ki]) for Ki in Kscalar])
"""


class Transport(object):
    def __init__(
        self,
        rhow,
        b,
        rhocp,
        muf,
        porosity,
        perm,
        conduct,
        onecomp,
        exact_sol,
        grid,
        pleft,
        pright,
        a,
        CL,
        CR,
    ):

        self.simulation = ComPASS.load_eos("linear_water")

        fluid_properties = self.simulation.get_fluid_properties()
        fluid_properties.specific_mass = rhow
        fluid_properties.volumetric_heat_capacity = b
        fluid_properties.dynamic_viscosity = muf
        self.simulation.set_rock_volumetric_heat_capacity(rhocp)

        self.on_the_left = lambda x: x <= grid.origin[0]
        self.on_the_left_in_the_middle = (
            lambda x, y: (x <= grid.origin[0])
            & (grid.extent[1] / 2 - a < y)
            & (y < grid.extent[1] / 2 + a)
        )
        self.on_the_right = lambda x: x >= grid.origin[0] + grid.extent[0]

        # only sequential allowed here
        assert ComPASS.mpi.master_proc_rank == ComPASS.mpi.proc_rank
        ComPASS.set_output_directory_and_logfile(__file__)

        self.simulation.init(
            grid=grid,
            set_dirichlet_nodes=self.select_dirichlet_nodes,
            cell_porosity=porosity,
            cell_permeability=perm,
            cell_thermal_conductivity=self.cell_diffusion_tensor_factory(
                conduct
            ),  # variable_conductivity, #conduct,
        )

        self.nb_nodes = self.simulation.global_number_of_nodes()
        self.nb_cells = self.simulation.global_number_of_cells()
        self.nb_points = self.nb_nodes + self.nb_cells
        self.CL = CL
        self.CR = CR
        self.grid = grid

        self.set_initial_values(pright, onecomp)
        self.set_boundary_conditions(pleft, pright, onecomp)

        #%%-- The three following functions are here just to hide
        #     the fact that we use temperatures as concentrations

    def retrieve_concentrations(self):
        # copy needed
        return np.copy(self.simulation.all_states().T)

    def set_concentrations(self, C):
        self.simulation.all_states().T[:] = C

    def set_source_term(self, S):
        # WARNING: only porous volume at dof not including rock volume
        porous_volume = self.simulation.all_Fourier_porous_volumes()
        thermal_sources = self.simulation.all_thermal_sources()
        thermal_sources[:] = porous_volume * S

    def cell_diffusion_tensor_factory(self, conduct):
        def cell_diffusion_tensor():
            cell_centers = self.simulation.compute_global_cell_centers()
            xc = cell_centers[:, 0]
            yc = cell_centers[:, 1]
            nbcells = cell_centers.shape[0]

            # tensor array
            celldiffusion = np.empty((nbcells, 3, 3), dtype=np.double)
            alpha_L = conduct
            alpha_T = conduct

            normVelocity = 0.277778e-05
            zeroA = np.zeros(nbcells)
            onesA = np.ones(nbcells)

            Dxx = alpha_L * normVelocity * onesA
            Dxy = zeroA
            Dyx = Dxy
            Dyy = alpha_T * normVelocity * onesA

            diff = np.array(
                [[Dxx, Dxy, zeroA], [Dyx, Dyy, zeroA], [zeroA, zeroA, zeroA]],
                dtype=np.double,
            )
            return diff.T

        return cell_diffusion_tensor

    def select_dirichlet_nodes(self):
        x = self.simulation.global_vertices()[:, 0]
        return self.on_the_left(x) | self.on_the_right(x)

    def set_boundary_conditions(self, pleft, pright, onecomp):
        def set_states(states, x, y):
            left = self.on_the_left(x)
            right = self.on_the_right(x)
            leftmid = self.on_the_left_in_the_middle(x, y)
            both = left | right

            states.p[left] = pleft
            states.p[right] = pright

            states.context[both] = 1
            states.S[both] = 1
            if onecomp:
                states.C[both] = 1.0
            else:
                states.C[left] = (1, 0)  # (0., 1)
                states.C[leftmid] = (0.0, 1.0)
                states.C[right] = (1, 0)

        verts = self.simulation.vertices()
        set_states(self.simulation.dirichlet_node_states(), verts[:, 0], verts[:, 1])

    def set_states_inj(self, states, x, y, idx):
        left = self.on_the_left(x)
        right = self.on_the_right(x)
        leftmid = self.on_the_left_in_the_middle(x, y)
        both = left | right

        states.T[right] = self.CR[idx]
        states.T[left] = self.CR[idx]
        states.T[leftmid] = self.CL[idx]

    def set_initial_values(self, p0, onecomp):
        def set_states(states, x):
            states.context[:] = 1
            states.p[:] = p0  # pleft + (pright - pleft) * (x - grid.origin[0]) / Lx
            states.T[:] = 0  # Tright
            states.S[:] = 1
            if onecomp:
                states.C[:] = 1.0
            else:
                states.C[:] = (1, 0)

        verts = self.simulation.vertices()
        cellcenters = self.simulation.compute_cell_centers()
        set_states(self.simulation.node_states(), verts[:, 0])
        set_states(self.simulation.cell_states(), cellcenters[:, 0])

    def plot_pressure(self):

        xy = self.simulation.compute_cell_centers()[:, 0:2]
        nx = self.grid.shape[0]
        ny = self.grid.shape[1]
        XX = xy[:, 0].reshape(ny, nx)
        YY = xy[:, 1].reshape(ny, nx)
        fig = plt.figure(1)
        cs0 = plt.contourf(
            XX,
            YY,
            np.reshape(self.simulation.cell_states().p, [ny, nx]),
            15,
            cmap="jet",
        )
        fig.colorbar(cs0)
        plt.title("Pressure")
        plt.xlabel("x in meters")
        plt.ylabel("y in meters")
        plt.draw()
        plt.savefig(ComPASS.to_output_directory("Pressure"), format="png")
        plt.pause(5)
        plt.clf()

    def plot_3d_concentrations(self, t, conc):

        Caq = ["Ca2+", "Na+", "Cl-", "K+"]
        xyz = self.simulation.compute_cell_centers()[:, 0:3]
        nx = self.grid.shape[0]
        ny = self.grid.shape[1]
        nz = self.grid.shape[2]
        XX = xyz[:, 0].reshape(nz, ny, nx)
        YY = xyz[:, 1].reshape(nz, ny, nx)
        ZZ = xyz[:, 2].reshape(nz, ny, nx)
        fig = plt.figure(1)
        fig.clf()

        for i in range(len(Caq)):
            Ci_nodes = conc[i][0 : self.nb_nodes]
            Ci_cells = conc[i][self.nb_nodes : self.nb_points]

            ax = plt.subplot(4, 1, i + 1)
            # lvs = [1e-10, 0.000147, 0.000294, 0.000442, 0.000589]
            # cs =plt.contourf(XX,YY,np.reshape(Ci_cells,[ny,nx]), 15, levels=lvs, vmin=1e-10, vmax=0.000589,cmap='jet')
            # cs = plt.contourf(XX, YY, ZZ, np.reshape(Ci_cells, [nz, ny, nx]), 15, cmap="jet")

            print("XX.size = ", XX.size)
            cs = mlab.contour3d(XX, YY, ZZ, np.reshape(Ci_cells, [nz, ny, nx]))
            # fig.colorbar(cs)
            plt.setp(ax.get_xticklabels(), visible=False)
            plt.gcf().subplots_adjust(hspace=0.4)
            plt.title(Caq[i])
            plt.suptitle("t=" + str(t / 3600.0) + "  (hrs)")
        plt.setp(ax.get_xticklabels(), visible=True)
        plt.xlabel("x in meters")
        plt.ylabel("y in meters")
        # plt.zlabel("z in meters")
        plt.draw()
        plt.pause(0.1)
        # if ((t==0)|(t==21600)|(t==28800)) :
        #     plt.savefig(ComPASS.to_output_directory('diff_cell_conc_at_time_'+str(t)),format='png')
        plt.savefig(
            ComPASS.to_output_directory("cell_conc_at_time_" + str(t)), format="png"
        )

    def plot_concentrations(self, t, conc):

        Caq = ["Ca2+", "Na+", "Cl-", "K+"]
        xy = self.simulation.compute_cell_centers()[:, 0:2]
        nx = self.grid.shape[0]
        ny = self.grid.shape[1]
        XX = xy[:, 0].reshape(ny, nx)
        YY = xy[:, 1].reshape(ny, nx)
        fig = plt.figure(1)
        fig.clf()

        for i in range(len(Caq)):
            Ci_nodes = conc[i][0 : self.nb_nodes]
            Ci_cells = conc[i][self.nb_nodes : self.nb_points]

            ax = plt.subplot(4, 1, i + 1)
            # lvs = [1e-10, 0.000147, 0.000294, 0.000442, 0.000589]
            # cs =plt.contourf(XX,YY,np.reshape(Ci_cells,[ny,nx]), 15, levels=lvs, vmin=1e-10, vmax=0.000589,cmap='jet')
            cs = plt.contourf(XX, YY, np.reshape(Ci_cells, [ny, nx]), 15, cmap="jet")
            fig.colorbar(cs)
            plt.setp(ax.get_xticklabels(), visible=False)
            plt.gcf().subplots_adjust(hspace=0.4)
            plt.title(Caq[i])
            plt.suptitle("t=" + str(t / 3600.0) + "  (hrs)")

            # Show the major grid lines with dark grey lines
            # plt.grid(b=True, which='major', color='#666666', linestyle='-')
            # Show the minor grid lines with very faint and almost transparent grey lines
            # plt.minorticks_on()
            # cs =plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)

        plt.setp(ax.get_xticklabels(), visible=True)
        plt.xlabel("x in meters")
        plt.ylabel("y in meters")
        plt.draw()
        # plt.pause(0.1)
        # if ((t==0)|(t==21600)|(t==28800)) :
        # 	plt.savefig(ComPASS.to_output_directory('diff_cell_conc_at_time_'+str(t)),format='png')
        plt.savefig(ComPASS.to_output_directory("cell_conc_at_time_" + str(t) + ".png"))

    def plot_1D_concentrations(self, t, conc):

        Caq = ["Ca2+", "Na+", "Cl-", "K+"]
        x = self.simulation.cell_centers()[:, 0]
        nx = self.grid.shape[0]
        ny = self.grid.shape[1]
        ligne = nx * int(ny / 2)
        fig2 = plt.figure(2)
        fig2.clf()

        for i in range(len(Caq)):
            Ci_nodes = conc[i][0 : self.nb_nodes]
            Ci_cells = conc[i][self.nb_nodes : self.nb_points]

            ax = plt.subplot(2, 1, 1)
            # cs =plt.plot(x[ligne:ligne+30], Ci_cells[ligne:ligne+30], label=Caq[i])
            cs = plt.plot(
                x[ligne : ligne + nx],
                Ci_nodes[ligne + int(ny / 2) : ligne + nx + int(ny / 2)],
                label=Caq[i],
            )
            plt.title("plot at nodes")
            plt.suptitle("t=" + str(t / 3600.0) + " (hrs)")
            plt.setp(ax.get_xticklabels(), visible=False)

            ax2 = plt.subplot(2, 1, 2)
            cs = plt.plot(
                x[ligne : ligne + nx],
                Ci_cells[ligne + int(ny / 2) : ligne + nx + int(ny / 2)],
                label=Caq[i],
            )
            plt.title("plot at cells")
            # plt.title("t=" + str(t / 3600.0) + " (hrs)")
            plt.setp(ax.get_xticklabels(), visible=True)
            # Show the major grid lines with dark grey lines
            # plt.grid(b=True, which='major', color='#666666', linestyle='-')
            # Show the minor grid lines with very faint and almost transparent grey lines
            # plt.minorticks_on()
            # cs =plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)

        plt.xlabel("x in meters")
        plt.ylabel("conc")
        plt.legend()
        plt.draw()
        # plt.pause(0.1)
        plt.savefig(
            ComPASS.to_output_directory("1D_cell_conc_at_time_" + str(t) + ".png")
        )
        # plt.savefig(ComPASS.to_output_directory('1D_cell_conc_at_time_%03d.png'%t))
        # if ((t%500)<1e-6 ) :
        # plt.savefig(ComPASS.to_output_directory('cell_conc_at_time_'+str(math.floor(t))),format='png')
        #    plt.savefig(ComPASS.to_output_directory('1D_cell_conc_at_time_'+str(t)+'.png'))

        # if ((t==0)|(t==21600)|(t==28800)) :
        # 	plt.savefig(ComPASS.to_output_directory('1D_cell_conc_at_time_'+str(t)),format='png')

    # the function that advects and diffuses all concentration
    def transport_concentrations(self, t, ts_manager, Cold, srcF):
        Nc = Cold.shape[0]
        Cnew = np.zeros_like(Cold)

        compute_all_concentrations = True

        newton = Newton(
            self.simulation,
            1e-5,
            30,
            linear_solver(simulation, legacy=True, tolerance=1e-6, max_iterations=150),
        )
        context = SimulationContext()
        # self.plot_pressure()

        # while loop to find suitable dt (other strategies are possible...)
        while compute_all_concentrations:

            for i in range(Nc):
                print("espece ========> ", i)
                # print ("ts_manager.current_step ========> ", ts_manager.current_step)

                self.set_source_term(srcF[i, :] / ts_manager.current_step)

                self.set_concentrations(Cold[i, :])
                self.set_states_inj(
                    self.simulation.dirichlet_node_states(),
                    self.simulation.vertices()[:, 0],
                    self.simulation.vertices()[:, 1],
                    i,
                )

                deltat = timestep.make_one_timestep(
                    newton, ts_manager.steps(), simulation_context=context
                )

                Cnew[i, :] = self.retrieve_concentrations()

            else:
                compute_all_concentrations = (
                    False  # job is done all concentrations have been computed with dt
                )
        # self.plot_1D_concentrations(t, Cnew)

        return Cnew
