import matplotlib.pyplot as plt
import numpy as np
from numpy import linalg as lg
from scipy.linalg import lstsq
from scipy import optimize
import math as ma

class Chemistry (object) :
	def __init__(self, T, W, c0, s0, S, A, B,  logKx, logKy, impflag,fictflag):

        	# Arguments:
        	#    T, W    : Total analytic concentrations for mobile (resp immobile)
        	#              components
        	#    c0, s0  : initial iterates for mobile (resp. immobile) components
        	#    S, A, B : Stoechiometric matrices (see Erhel/Kedrn report for
        	#              notations)
        	#    Kx, Ky  : vectors of equilibrium constants for aqueous (resp ion
        	#              exchange) reactions
        	#    impflag : array of flags for imposed concentration (impflag(j)=0 
        	#              means concentration of component j is known. In that case,
        	#              the total concentration must contain the given value
        	#    fictflag: array of flags to take into account ion exchange or not...  
          
		self.rhs = np.append(T, W)
		if (A.shape[0] != B.shape[0]):
			print('Wrong size for A or B');
		self.StMat = np.block([[S, np.zeros((S.shape[0], B.shape[1]))], [A, B]])  
		self.logK =np.append(logKx, logKy)
		self.impflag = impflag
		self.fictflag = fictflag
		
		for i in range(len(T)) :
			if c0[i]==0.0 :
				c0[i] = 1e-8
		for i in range(len(W)) :
			if s0[i]==0.0 :
				s0[i] = 1e-8

		self.zin = np.append(np.log(c0), np.log(s0))
		
		# useful numbers for post-process mobile and fixed concentrations
		self._Nc = len(T)
		self._Ns = len(W)
		# added for for post-process
		self.S = S
		self.A = A
		self.B = B
		self.logKx = logKx
		self.logKy = logKy
		# added
		self.T = T
		self.W = W

	def fixed_free (self):

	        # package things before calling solver
		fixed= np.flatnonzero(self.impflag==0);   
		free = np.flatnonzero(self.impflag!=0);  

	        # eliminate fixed concentrations
		if (len(fixed)>0):
			logKff =  np.transpose([self.logK]) + self.StMat[:, fixed]*np.log(self.rhs[fixed]);
		else :  logKff = np.transpose([self.logK]) 

		self.fictflag=self.fictflag[free];
		self.zin=self.zin[free];
		        
		return fixed, free, logKff
	def get_Nc (self):
   		return self._Nc
   
	def fJacChimie(self, z, rhs, StMat, logKff):

		z=z.reshape(z.size,1)
		rhs=rhs.reshape(rhs.size,1)
		fictflag = self.fictflag.reshape(self.fictflag.size,1)
		
		big=50; small=-500;
  
		tmp=np.exp(np.minimum(big, np.maximum(small, logKff+np.matmul(StMat, z))))

		cs = np.exp(np.minimum(big, np.maximum(small,z)));
	
		fcs = fictflag*cs
		fcs=fcs.reshape(1, fcs.size)
		f = rhs - fictflag * cs - np.matmul(np.transpose(StMat), tmp); 
   
		jac = - np.diag(fcs[0]) - np.matmul(np.matmul(np.transpose(StMat), np.diag(tmp.reshape(1, tmp.size)[0])), StMat);
	
		f=f.reshape(1, f.size)

		return f[0], jac

	def  solve (self) :    
        	# SOLVEQ - solves chemical equilibrium systems 
        	#          for aquatic chemistry
        	#          Uses scipy.optimize.root
     
		[fixed, free, logKff] = self.fixed_free ()
		
		zsol = optimize.root(self.fJacChimie, self.zin, args = (self.rhs[free], self.StMat[:, free], logKff), jac=True, method='lm', options={'ftol': 1.49012e-08});
		#zsol = optimize.root(self.fJacChimie, self.zin,  args = (self.rhs[free], self.StMat[:, free], logKff), jac=None, method='krylov', options={'disp': True})  #broyden2
            
		#print ('=====> zsol Result = ', zsol)
		#print ('%%%%%%%%%%%%%> zsol.success = ', zsol.success)
		#print ('=====> zsol.x = ', zsol.x);
		#print ('==> zsol.status = ', zsol.status) 
		#print ('=====> zsol.nfev = ', zsol.nfev)
		#print ('==> zsol.njev = ', zsol.njev)
		#print ('=====> zsol.nit = ', zsol.nit)
	        #print ('=====> zsol.fun = ', zsol.fun)
        
		status = zsol.success
	        #nfuneval = zsol.nfev 
		niter = zsol.njev  #zsol.nit number of iterations #TO DO : to get from the program as result specially for 'lm' and 'hybr' methods 
                  #because there is no nit attribute for optimizeResult objet not like other methods

		zzsol=np.zeros_like(self.impflag, dtype=float);     
		
		zzsol[free]=zsol.x
		zzsol[fixed]=np.log(self.rhs[fixed])
						         
		return zzsol, status, niter 

	def  post_process (self, zzsol) :
		
	        # post-process to get mobile and fixed concentrations
		Nc = self._Nc; Ns = self._Ns
		zzsol=zzsol.reshape(zzsol.size,1)

		logcsol=zzsol[0:Nc];
		logssol=zzsol[Nc:Nc+Ns];
		#print("logcsol=", logcsol)
		csol=np.exp(logcsol);
		ssol=np.exp(logssol);

		fixed= np.flatnonzero(self.impflag==0); 
		logKx = self.logKx.reshape(self.logKx.size,1)
		logKy = self.logKy.reshape(self.logKy.size,1)
		
		logx = np.matmul(self.S,logcsol)+logKx; 
		x=np.exp(logx);
		logy = np.matmul(self.A,logcsol)+np.matmul(self.B,logssol)+logKy; 
		y=np.exp(logy);
		C = csol+np.matmul(np.transpose(self.S), x)  
		C[fixed] = csol[fixed]; 

		F = np.matmul(np.transpose(self.A), y) 
		F[fixed]=0;

		csol=csol.reshape(1, csol.size)[0]
		ssol=ssol.reshape(1, ssol.size)[0]
		
		##### calcul de psi_prime_dT ########################################
		[fixed, free, logKff] = self.fixed_free ()
		[ff, jac]=self.fJacChimie(zzsol[free], self.rhs[free], self.StMat[:, free], logKff)
		sbr=np.block([[np.eye(Nc, dtype=float)], [np.zeros((1, Nc))]])
		#D=np.diag(1./np.sum(abs(jac),1));
		#X= np.linalg.solve(-np.matmul(D,jac), np.matmul(D,sbr)) 
		
		#X = np.linalg.solve(-jac, sbr)
		X, res, rk, s =  lstsq(-jac, sbr)  #### scipy.linalg.lstsq
		psi_prime_dT=np.matmul(np.matmul(np.matmul(np.transpose(self.A), np.diagflat(y)),np.block([[self.A,self.B]])), X)

		return C, F, csol, ssol, psi_prime_dT






