#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 11:36:10 2019

@author: kern
"""

import ComPASS
from ComPASS.utils.units import *
from ComPASS.timeloops import standard_loop
from ComPASS import newton

import numpy as np


def CrudeNewton():
    return newton.Newton(1e-6, 8, newton.LinearSolver(1e-6, 150))


Lx, Ly = 2.1, 1.0
nx, ny = 100, 50

pout = 1
fluxInj1 = 2.25e-2
fluxInj2 = 2.25e-2
Tinj1 = 0.3
Tinj2 = 0.3
omegaA, omegaB = 0.25, 0.5
permA, permB = 1e-2, 1e-5
diffusion = 2e-3

rhow = 1  # 1E3
b = 1  # 4.2e+3          # specific heat in J/kg/K
muf = 1  # 3E-4

ComPASS.load_eos("linear_water")
fluid_properties = ComPASS.get_fluid_properties()
fluid_properties.specific_mass = rhow
fluid_properties.volumetric_heat_capacity = b
fluid_properties.dynamic_viscosity = muf
ComPASS.set_rock_volumetric_heat_capacity(0)

ComPASS.set_gravity(0)

grid = ComPASS.Grid(shape=(nx, ny, 1), extent=(Lx, Ly, Lx / nx), origin=(0, 0, 0))

inflow1 = (
    lambda x, y: (y >= grid.origin[1] + grid.extent[1])
    & (grid.origin[0] + 0.3 <= x)
    & (x <= grid.origin[0] + 0.4)
)
inflow2 = (
    lambda x, y: (x <= grid.origin[0])
    & (grid.origin[1] + 0.3 <= y)
    & (y <= grid.origin[1] + 0.4)
)
outflow = (
    lambda x, y: (x >= grid.origin[0] + grid.extent[0])
    & (grid.origin[1] + 0.2 <= y)
    & (y <= grid.origin[1] + 0.4)
)


def select_dirichlet_Pnodes():
    x = ComPASS.global_vertices()[:, 0]
    y = ComPASS.global_vertices()[:, 1]
    return outflow(x, y)


def select_dirichlet_Tnodes():
    x = ComPASS.global_vertices()[:, 0]
    y = ComPASS.global_vertices()[:, 1]
    return inflow1(x, y) | inflow2(x, y)


def cell_permeability_factory(grid):
    def cell_permeability():
        cell_centers = ComPASS.compute_global_cell_centers()
        xc = cell_centers[:, 0]
        yc = cell_centers[:, 1]
        nbcells = cell_centers.shape[0]
        obstacle = (
            (grid.origin[0] + 1 <= xc)
            & (xc <= grid.origin[0] + 1.1)
            & (yc <= grid.origin[1] + 0.9)
        )
        cellperm = permA * np.ones((nbcells), dtype=np.double)
        cellperm[obstacle] = permB
        return cellperm

    return cell_permeability


def cell_porosity_factory(grid):
    def cell_porosity():
        cell_centers = ComPASS.compute_global_cell_centers()
        xc = cell_centers[:, 0]
        yc = cell_centers[:, 1]
        nbcells = cell_centers.shape[0]
        obstacle = (
            (grid.origin[0] + 1 <= xc)
            & (xc <= grid.origin[0] + 1.1)
            & (yc <= grid.origin[1] + 0.9)
        )
        cellporos = omegaA * np.ones((nbcells), dtype=np.double)
        cellporos[obstacle] = omegaB
        return cellporos

    return cell_porosity


ComPASS.set_output_directory_and_logfile(__file__)

ComPASS.init(
    mesh=grid,
    cell_porosity=cell_porosity_factory(grid),
    cell_permeability=cell_permeability_factory(grid),
    set_pressure_dirichlet_nodes=select_dirichlet_Pnodes,
    set_temperature_dirichlet_nodes=select_dirichlet_Tnodes,
    cell_thermal_conductivity=diffusion,
)


def set_boundary_conditions():
    def set_states(states, x, y):
        in1 = inflow1(x, y)
        states.T[in1] = Tinj1
        in2 = inflow2(x, y)
        states.T[in2] = Tinj2
        out = outflow(x, y)
        states.p[out] = pout
        allbc = in1 | in2 | out
        states.context[allbc] = 1
        states.S[allbc] = 1
        states.C[allbc] = 1.0

    verts = ComPASS.vertices()
    set_states(ComPASS.dirichlet_node_states(), verts[:, 0], verts[:, 1])


def set_initial_values():
    def set_states(states, x):
        states.context[:] = 1
        states.p[:] = pout  # pleft + (pright - pleft) * (x - grid.origin[0]) / Lx
        states.T[:] = 0
        states.S[:] = 1
        states.C[:] = 1.0

    verts = ComPASS.vertices()
    cellcenters = ComPASS.compute_cell_centers()
    set_states(ComPASS.node_states(), verts[:, 0])
    set_states(ComPASS.cell_states(), cellcenters[:, 0])


set_initial_values()
set_boundary_conditions()


def set_boundary_flux():
    face_centers = ComPASS.face_centers()
    x = face_centers[:, 0]
    y = face_centers[:, 1]
    # injection 1 on top
    Neumann = ComPASS.NeumannBC()
    Neumann.molar_flux[:] = fluxInj1
    ComPASS.set_Neumann_faces(inflow1(x, y), Neumann)
    # injection 2 on left side
    Neumann = ComPASS.NeumannBC()
    Neumann.molar_flux[:] = fluxInj2
    ComPASS.set_Neumann_faces(inflow2(x, y), Neumann)


set_boundary_flux()


final_time = 60
nitermax = 100
dt = 1e-2
standard_loop(final_time=final_time, output_period=2, initial_timestep=dt)
# standard_loop(nitermax=nitermax, output_every=1, output_period=10*dt, initial_timestep=dt, newton=CrudeNewton())
