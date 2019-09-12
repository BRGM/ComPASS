import ComPASS

# import doublet_utils
from ComPASS.utils.units import *
from ComPASS.timeloops import standard_loop
import matplotlib.pyplot as plt
import numpy as np

# from numpy import linalg as lg
from scipy.linalg import lstsq
from scipy import optimize, sparse
import transport as tr
import chemistry as ch
from algo import Struct
import math as ma

# import Newton_krylov as NK
from NonLinSolve import nonlin_solve
from scipy.sparse.linalg import LinearOperator


#######################################################################################
class formulation(Struct):
    nom = str


####################   class evoltrch  ################################################
class Coupling(object):
    def __init__(self, Chobj, TrObj, nb_time_steps):

        self.Tot = np.zeros((Chobj._Nc, TrObj.nb_points))
        self.W = np.zeros((Chobj._Ns, TrObj.nb_points))
        self.C = np.zeros_like(self.Tot)
        self.F = np.zeros_like(self.Tot)
        self.Csol = np.zeros((nb_time_steps, Chobj._Nc))

        # initialisations
        for i in range(Chobj._Nc):
            self.Tot[i][:] = Chobj.T[i]
        for i in range(Chobj._Ns):
            self.W[i][:] = Chobj.W[i]

        # separation to mobile and immobile constrations
        for ix in range(TrObj.nb_points):

            Chobj.rhs = np.append(self.Tot[:, ix], self.W[:, ix])
            [sol, status, niter] = Chobj.solve()
            [Caux, Faux, cc, ss, psi] = Chobj.post_process(sol)

            self.C[:, ix] = Caux.T
            self.F[:, ix] = Faux.T

        self.jacCh = np.zeros(
            ((TrObj.nb_points) * Chobj._Nc, (TrObj.nb_points) * Chobj._Nc)
        )

    ############### JacVec   ######################################################
    """   #w = v - jac_ch*v+ jac_ch*(jac_tr\(diag(phi)*v));
	   
	def JacVec(self, v, jac_ch, jac_tr):

	     sm = np.matmul(np.diagflat(tr.omega_reservoir*np.ones((tr.nb_points)*Nc)), v)
	     #D=np.diag(1./np.sum(abs(jact),1))
	     #tmp= np.linalg.solve(np.matmul(D,jact), np.matmul(D,sm)) 
	     
	     #tmp = np.linalg.solve(jact, sm)
	     tmp = np.matmul(np.linalg.pinv(jac_tr), sm)

	     r= v - np.matmul(jac_ch, v) + np.matmul(jac_ch, tmp)
	     
	     return r
	"""
    ### F formulation for coupling transport with chemistry #######################
    def meth_F(self, Y, Cold, Fold, csol, Chobj, Trobj, t, ts_manager):
        nb_points = Trobj.nb_points
        Nc = Chobj._Nc
        F = Y.reshape(Chobj._Nc, nb_points)

        ff = np.zeros_like(Fold)
        g = np.zeros_like(Fold)

        # transport pour chaque composante non imposee
        srcF = Fold - F
        # print("srcF = ", srcF)
        Cnew = Trobj.transport_concentrations(t, ts_manager, Cold, srcF)

        self.C = Cnew
        # Trobj.plot_concentrations(t, self.C)
        g = Cnew + F

        # chimie pour chaque point du maillage
        """
		for ix in range(nb_points) :
		    #Chobj.c0 = [1e-6, 1e-6, 1e-6, 1e-6] 
		    #Chobj.s0 = [0.54849]
		    Chobj.c0 = csol[:, ix-1]
		    Chobj.rhs = np.append(g[:, ix], self.W[:, ix])
		    [sol, status, niter] = Chobj.solve()
		    [Caux, Faux, cc, ss, psi_prime_dT]=Chobj.post_process (sol)
		    ff[:,ix]=-Faux.T + F[:,ix]
		    csol[:,ix]=cc.T
		    self.jacCh[ix::nb_points,ix::nb_points] = psi_prime_dT
		"""

        Chobj.rhs = np.append(g[:, 0], self.W[:, 0])
        [sol, status, niter] = Chobj.solve()
        [Caux, Faux, cc, ss, psi_prime_dT] = Chobj.post_process(sol)
        ff[:, 0] = -Faux.T + F[:, 0]
        csol[:, 0] = cc.T
        self.jacCh[0::nb_points, 0::nb_points] = psi_prime_dT

        for ix in range(nb_points):
            # Chobj.c0 = [1e-6, 1e-6, 1e-6, 1e-6]
            # Chobj.s0 = [0.54849]

            cpt = 0
            if ix > 0:
                Chobj.c0 = csol[:, ix - 1]
                for ic in range(Nc):
                    if abs(g[ic, ix] - g[ic, ix - 1]) < 1e-6:
                        cpt = cpt + 1

                if cpt == Nc:
                    ff[:, ix] = ff[:, ix - 1]
                    self.jacCh[ix::nb_points, ix::nb_points] = self.jacCh[
                        ix - 1 :: nb_points, ix - 1 :: nb_points
                    ]
                    csol[:, ix] = csol[:, ix - 1]
                else:
                    Chobj.rhs = np.append(g[:, ix], self.W[:, ix])
                    [sol, status, niter] = Chobj.solve()
                    [Caux, Faux, cc, ss, psi_prime_dT] = Chobj.post_process(sol)
                    ff[:, ix] = -Faux.T + F[:, ix]
                    csol[:, ix] = cc.T
                    self.jacCh[ix::nb_points, ix::nb_points] = psi_prime_dT

        return (
            ff.flatten(),
            Cnew.flatten(),
            self.jacCh,
        )  ### pas besoin de retourner jacCh si on le garde comme attribut de classe

    #######  CF formulation for coupling transport with chemistry  ########################################
    def meth_CF(self, Y, Cold, Fold, csol, Chobj, Trobj, ts_manager):
        nb_points = Trobj.nb_points
        if 2 * Chobj._Nc * nb_points != len(Y):
            print("wrong size for Y")

        C = np.zeros_like(Fold)
        F = np.zeros_like(Fold)
        f1 = np.zeros_like(Fold)
        f2 = np.zeros_like(Fold)

        tmp = Y.reshape(2 * Chobj._Nc, nb_points)
        for i in range(Chobj._Nc):
            C[i, :] = tmp[i][0:nb_points]
            F[i, :] = tmp[i + Chobj._Nc][0:nb_points]

        tot = C + F

        # transport pour chaque composante non imposee
        srcF = Fold - F
        Cnew = Trobj.transport_concentrations(ts_manager, Cold, srcF)
        ###self.C = Cnew ?? faut il mettre a jour self.C ici ?? a verifier

        f1 = Cnew

        # chimie pour chaque point du maillage
        """
	 	Chobj.rhs = np.append(tot[:, 1], self.W[:, 1])  
	 	[sol, status, niter] = Chobj.solve()
	 	[Caux, Faux, cc, ss, psi_prime_dT]=Chobj.post_process (sol)
	 	f2[:, 1]= F[:, 1] -Faux.T 
	 	csol[:,1]=cc.T
	 	self.jacCh[1::nb_points,1::nb_points] = psi_prime_dT
	 	"""
        for ix in range(nb_points):
            # Chobj.c0 = csol[:, ix]
            Chobj.c0 = [1e-6, 1e-6, 1e-6, 1e-6]
            Chobj.s0 = [0.54849]
            Chobj.rhs = np.append(tot[:, ix], self.W[:, ix])
            [sol, status, niter] = Chobj.solve()
            [Caux, Faux, cc, ss, psi_prime_dT] = Chobj.post_process(sol)
            f2[:, ix] = F[:, ix] - Faux.T
            csol[:, ix] = cc.T
            self.jacCh[ix::nb_points, ix::nb_points] = psi_prime_dT
            """		    
		    
		    if ix>1 :
		    	#Chobj.c0 = csol[:, ix-1]
		    	if (tot[:,ix].all() == tot[:,ix-1].all()) :
		    		f2[:,ix] = f2[:,ix-1]
		    		self.jacCh[ix::nb_points,ix::nb_points] = self.jacCh[ix-1::nb_points,ix-1::nb_points]
		    		csol[:,ix]=csol[:,ix-1]
		    	else :
		    		Chobj.rhs = np.append(tot[:, ix], self.W[:, ix])  
		    		[sol, status, niter] = Chobj.solve()
		    		[Caux, Faux, cc, ss, psi_prime_dT]=Chobj.post_process (sol)
		    		f2[:,ix]= F[:,ix] -Faux.T 
		    		csol[:,ix]=cc.T
		    		self.jacCh[ix::nb_points,ix::nb_points] = psi_prime_dT
	 		"""
        Ysor = np.append(f1.flatten(), f2.flatten())

        return Ysor, f1.flatten(), self.jacCh

    ###  CFT formulation for coupling transport with chemistry  ################################################
    def meth_CFT(self, Y, Cold, Fold, csol, Chobj, Trobj, ts_manager):

        Nc = Chobj._Nc
        nb_points = Trobj.nb_points
        if 3 * Nc * nb_points != len(Y):
            print("wrong size for Y")

        C = np.zeros_like(Fold)
        Tot = np.zeros_like(Fold)
        F = np.zeros_like(Fold)
        f1 = np.zeros_like(Fold)
        f2 = np.zeros_like(Fold)
        f3 = np.zeros_like(Fold)

        tmp = Y.reshape(3 * Nc, nb_points)
        for i in range(Nc):
            C[i, :] = tmp[i][0:nb_points]
            Tot[i, :] = tmp[i + Nc][0:nb_points]
            F[i, :] = tmp[i + 2 * Nc][0:nb_points]

        f2 = Tot - C - F

        # transport pour chaque composante non imposee
        srcF = Fold - F
        [Cnew, dt] = Trobj.transport_concentrations(ts_manager, Cold, srcF)

        f1 = Cnew

        # chimie pour chaque point du maillage
        for ix in range(nb_points):
            # Chobj.c0 = csol[:, ix]
            Chobj.c0 = [1e-6, 1e-6, 1e-6, 1e-6]
            Chobj.s0 = [0.54849]

            Chobj.rhs = np.append(Tot[:, ix], self.W[:, ix])

            [sol, status, niter] = Chobj.solve()
            [Caux, Faux, cc, ss, psi_prime_dT] = Chobj.post_process(sol)

            f3[:, ix] = F[:, ix] - Faux.T
            csol[:, ix] = cc.T

            self.jacCh[ix::nb_points, ix::nb_points] = psi_prime_dT

        Ysor = np.append(f1.flatten(), np.append(f2.flatten(), f3.flatten()))

        return Ysor, f1.flatten(), self.jacCh

    ########### Simulation : evolution transport-chemistry (the time loop) ##############################
    def Simul_time_loop(self, xr, t, final_time, ts_manager, meth, Chobj, Trobj):

        # initial condition
        if meth.nom == "meth_CFT":
            Y0 = np.append(
                self.C.flatten(), np.append(self.Tot.flatten(), self.F.flatten())
            )
        elif meth.nom == "meth_CF":
            Y0 = np.append(self.C.flatten(), self.F.flatten())
        elif meth.nom == "meth_F":
            Y0 = self.F.flatten()
        else:
            error("Wrong value of formulation")

        Y = Y0
        csol = np.zeros_like(self.C)
        step = 0
        while t < final_time:  # boucle sur les pas de temps

            print("Doing time : ", t)

            Fold = self.F
            Cold = self.C

            if meth.nom == "meth_CFT":
                Y, info = nonlin_solve(
                    lambda Y: self.meth_CFT(
                        Y, Cold, Fold, csol, Chobj, Trobj, t, ts_manager
                    ),
                    Y,
                    jacobian="couplingJacCFT",
                    verbose=1,
                )
            elif meth.nom == "meth_CF":
                Y, info = nonlin_solve(
                    lambda Y: self.meth_CF(
                        Y, Cold, Fold, csol, Chobj, Trobj, t, ts_manager
                    ),
                    Y,
                    jacobian="couplingJacCF",
                    verbose=1,
                )
            elif meth.nom == "meth_F":
                Y, info = nonlin_solve(
                    lambda Y: self.meth_F(
                        Y, Cold, Fold, csol, Chobj, Trobj, t, ts_manager
                    ),
                    Y,
                    jacobian="couplingJacF",
                    verbose=1,
                )
                # Y = optimize.root(meth_F, Y, args = (Cfree,Ffree, t, dt, csol, ssol), jac=None, method='krylov');
                # Y = optimize.fixed_point(meth_F, Y, args = (Cfree,Ffree, t, dt, csol, ssol))
                # Y, info = NK.newton_krylov(lambda Y:self.meth_F(Y,Cold,Fold, t, dt,csol, Chobj, Trobj), Y, method='gmres', verbose=1)
            else:
                error("Wrong value of formulation")

            Nc = Chobj._Nc
            nb_points = Trobj.nb_points
            if meth.nom == "meth_CFT":
                tmp = Y.reshape(3 * Nc, nb_points)
                for i in range(Nc):
                    self.C[i, :] = tmp[i][0:nb_points]
                    self.Tot[i, :] = tmp[i + Nc][0:nb_points]
                    self.F[i, :] = tmp[i + 2 * Nc][0:nb_points]
            elif meth.nom == "meth_CF":
                tmp = Y.reshape(2 * Nc, nb_points)
                for i in range(Nc):
                    self.C[i, :] = tmp[i][0:nb_points]
                    self.F[i, :] = tmp[i + Nc][0:nb_points]
            elif meth.nom == "meth_F":
                self.F = Y.reshape(Nc, nb_points)
            else:
                error("Wrong value of formulation")

            Trobj.plot_1D_concentrations(t, self.C)
            Trobj.plot_concentrations(t, self.C)

            x = ComPASS.vertices()[:, 0]
            y = ComPASS.vertices()[:, 1]
            z = ComPASS.vertices()[:, 2]
            on_xr = (
                lambda x, y, z: (x == xr) & (y == Trobj.grid.extent[1] / 2) & (z == 0)
            )
            # on_xr = lambda x,y,z : (x == xr) & (y == 0) & (z == 0)   # test for ny =1

            self.Csol[step, :] = self.C[:, np.nonzero(on_xr(x, y, z))[0]].T

            t = t + ts_manager.current_step
            step = step + 1
