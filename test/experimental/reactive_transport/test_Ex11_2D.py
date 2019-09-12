import chemistry as ch
import transport as tr
import coupling as cp
import numpy as np
import ComPASS
import matplotlib.pyplot as plt
from ComPASS.timestep_management import FixedTimeStep, TimeStepManager

############################################################################
# exemple d'echange d'ions (example 11, doc phreeqc)
#
# 6 especes, aq=[[Ca2+, Na+, Cl-, K+], sol=[KX, CaX2, NaX]
#             x=[]; y=[CaX, NaX, KX]
# 4 composantes, c=[Ca2+, Na+, Cl-, K+], s=[X] <---  !!! composant fictif
############################################################################

# ======================= Chemical data ====================================#
S = np.zeros((0, 4))
A = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1]])
B = np.array([2, 1, 1])
B = B.reshape(B.size, 1)

logKx = np.zeros((0, 1))
# logKy = np.log(10)*np.array([0.8, 0, 0.7])
logKy = np.log(10) * np.array([3.4576, 0, 0.7])

T0 = [0.0, 1.5493e-3, 0.0, 7.507e-4]
# T0 = [1e-10, 1.5493e-3, 1e-10, 7.507e-4];
W0 = [1.1e-3]
impflag = np.array([1, 1, 1, 1, 1])
fictflag = np.array([1, 1, 1, 1, 0])
c0 = [1e-7, 1e-3, 1.0e-8, 1.0e-7]
s0 = [11.9e-3]
# ======================= Physical data ======================================#

rhow = 1
b = 1  # specific heat in J/kg/K
rhocp = 0  # volumetric heat capacity
muf = 1
porosity = 0.97  # 0.97 # reservoir porosity
perm = 1  # reservoir permeability in m^2
conduct = 5.55556e-9  # 0.002 diff     # bulk thermal conductivity in W/m/K
onecomp = True
exact_sol = True

pleft, pright = 0.277778e-05 * 0.12, 0  # 0.277778e-05*Lx, 0
# Tleft, Tright = degC2K(60), degC2K(100)
CL, CR = [6e-4, 0.0, 1.2e-3, 0.0], [0, 0.0, 0, 0.0]

# ======================= Geometry data =====================================#
a = 0.005  #
Lx, Ly, Lz = 0.12, 0.015, 0.003  # 0.12, 0.015, 0.004
nx, ny, nz = 60, 10, 1

# ===================== Chemical data recovery ==============================#
Chobj = ch.Chemistry(T0, W0, c0, s0, S, A, B, logKx, logKy, impflag, fictflag)

# ===================== Transport data recovery =============================#
grid = ComPASS.Grid(shape=(nx, ny, nz), extent=(Lx, Ly, Lz), origin=(0, 0, 0))
Trobj = tr.Transport(
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
)

# ===================== The simulation parameters  ==========================#
t = 0  # initial time
dt = 100  #  720 #tr.final_time/1e3
final_time = 72000  # 20 heures #29160 #

ts_manager = FixedTimeStep(dt)
# ts_manager = TimeStepManager(initial_timestep=720,)

tt = np.arange(t, final_time, dt)
nb_time_steps = tt.size
xr = 0.08

# choice of the method for coupling transport with chemistry
meth = cp.formulation(nom="meth_F")  #  meth_CF, meth_CFT

# ===================== Init coupled vectors T, C, F ========================#
evol = cp.Coupling(Chobj, Trobj, nb_time_steps)

# ===================== The simulation time loop  ===$=======================#
Csol = evol.Simul_time_loop(xr, t, final_time, ts_manager, meth, Chobj, Trobj)

# ===================== Plot of Elution curves  =============================#
# Trobj.plot_evol_concentrations(U*tt/xr, evol.Csol)

fig = plt.figure(3)
U = (pleft - pright) / Lx  # = 2.77778e-06

Caq = ["Ca2+", "Na+", "Cl-", "K+"]
for i in range(len(Caq)):
    plt.plot(U * tt / xr, evol.Csol[:, i], label=Caq[i])

plt.xlabel("pore volume")
plt.ylabel("concentrations")
plt.legend()
plt.draw()
plt.pause(0.1)
plt.savefig(ComPASS.to_output_directory("evol_conc_at_x=_0.08"), format="png")
plt.clf()

# plt.savefig(ComPASS.to_output_directory('evol_cell_conc_at_x='+str(xr)),format='png')
# ===================== End of simulation  =================================#
