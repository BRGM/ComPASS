
from __future__ import division, print_function, absolute_import

import sys
import numpy as np
from scipy._lib.six import callable, exec_, xrange
import  scipy.optimize.nonlin as nls
from scipy.linalg import norm, solve, inv, qr, svd, LinAlgError
from numpy import asarray, dot, vdot
import scipy.sparse.linalg
import scipy.sparse
from scipy.linalg import get_blas_funcs
import inspect
from scipy._lib._util import getargspec_no_self as _getargspec
#from .linesearch import scalar_search_wolfe1, scalar_search_armijo
from scipy.optimize.linesearch import scalar_search_wolfe1, scalar_search_armijo
import matplotlib.pyplot as plt
import ComPASS

__all__ = [
    'broyden1', 'broyden2', 'anderson', 'linearmixing',
    'diagbroyden', 'excitingmixing', 'newton_krylov']

#------------------------------------------------------------------------------
# Utility functions
#------------------------------------------------------------------------------


def nonlin_solve(F, x0, jacobian='krylov', iter=None, verbose=False,   
                 maxiter=None, f_tol=None, f_rtol=None, x_tol=None, x_rtol=None,
                 tol_norm=None, line_search='armijo', callback=None,
                 full_output=True, raise_exception=True):
    """
    Find a root of a function, in a way suitable for large-scale problems.

    Parameters
    ----------
    %(params_basic)s
    jacobian : Jacobian
        A Jacobian approximation: `Jacobian` object or something that
        `asjacobian` can transform to one. Alternatively, a string specifying
        which of the builtin Jacobian approximations to use:

            krylov, broyden1, broyden2, anderson
            diagbroyden, linearmixing, excitingmixing

    %(params_extra)s
    full_output : bool
        If true, returns a dictionary `info` containing convergence
        information.
    raise_exception : bool
        If True, a `NoConvergence` exception is raise if no solution is found.

    See Also
    --------
    asjacobian, Jacobian

    Notes
    -----
    This algorithm implements the inexact Newton method, with
    backtracking or full line searches. Several Jacobian
    approximations are available, including Krylov and Quasi-Newton
    methods.

    References
    ----------
    .. [KIM] C. T. Kelley, \"Iterative Methods for Linear and Nonlinear
       Equations\". Society for Industrial and Applied Mathematics. (1995)
       http://www.siam.org/books/kelley/fr16/index.php

    """

    condition = nls.TerminationCondition(f_tol=f_tol, f_rtol=f_rtol,
                                     x_tol=x_tol, x_rtol=x_rtol,
                                     iter=iter, norm=tol_norm)

    x0 = nls._as_inexact(x0)
    #func = lambda z: _as_inexact(F(_array_like(z, x0))).flatten()
    func = lambda z: F(nls._array_like(z, x0))
    x = x0.flatten()

    dx = np.inf
    Fx, _, _= func(x)
    #Fx = func(x)
    Fx_norm = norm(Fx)

    jacobian = asjacobian(jacobian)
    jacobian.setup(x.copy(), Fx, func)

    if maxiter is None:
        if iter is not None:
            maxiter = iter + 1
        else:
            maxiter = 100*(x.size+1)

    if line_search is True:
        line_search = 'armijo'
    elif line_search is False:
        line_search = None

    if line_search not in (None, 'armijo', 'wolfe'):
        raise ValueError("Invalid line search")

    # Solver tolerance selection
    gamma = 0.9
    eta_max = 0.9999
    eta_treshold = 0.1
    eta = 1e-3

    for n in xrange(maxiter):
        status = condition.check(Fx, x, dx)
        if status:
            break

        # The tolerance, as computed for scipy.sparse.linalg.* routines
        tol = min(eta, eta*Fx_norm)
        dx = -jacobian.solve(Fx, tol=tol)

        if norm(dx) == 0:
            raise ValueError("Jacobian inversion yielded zero vector. "
                             "This indicates a bug in the Jacobian "
                             "approximation.")

        # Line search, or Newton step
        if line_search:
            s, x, Fx, Fx_norm_new = _nonlin_line_search(func, x, Fx, dx,
                                                        line_search)
        else:
            s = 1.0
            x = x + dx
            Fx, _,_ = func(x)
            #Fx = func(x)
            Fx_norm_new = norm(Fx)

        jacobian.update(x.copy(), Fx)

        if callback:
            callback(x, Fx)

        # Adjust forcing parameters for inexact methods
        eta_A = gamma * Fx_norm_new**2 / Fx_norm**2
        if gamma * eta**2 < eta_treshold:
            eta = min(eta_max, eta_A)
        else:
            eta = min(eta_max, max(eta_A, gamma*eta**2))

        Fx_norm = Fx_norm_new

        # Print status
        if verbose:
            sys.stdout.write("%d:  |F(x)| = %g; step %g; tol %g\n" % (
                n, norm(Fx), s, eta))
            sys.stdout.flush()
    else:
        if raise_exception:
            raise NoConvergence(nls._array_like(x, x0))
        else:
            status = 2

    if full_output:
        info = {'nit': condition.iteration,
                'fun': Fx,
                'status': status,
                'success': status == 1,
                'message': {1: 'A solution was found at the specified '
                               'tolerance.',
                            2: 'The maximum number of iterations allowed '
                               'has been reached.'
                            }[status]
                }
        return nls._array_like(x, x0), info
    else:
        return nls._array_like(x, x0)

def _nonlin_line_search(func, x, Fx, dx, search_type='armijo', rdiff=1e-8,
                        smin=1e-2):
    tmp_s = [0]
    tmp_Fx = [Fx]
    tmp_phi = [norm(Fx)**2]
    s_norm = norm(x) / norm(dx)

    def phi(s, store=True):
        if s == tmp_s[0]:
            return tmp_phi[0]
        xt = x + s*dx
        v,_,_ = func(xt)
        p = nls._safe_norm(v)**2
        if store:
            tmp_s[0] = s
            tmp_phi[0] = p
            tmp_Fx[0] = v
        return p

    def derphi(s):
        ds = (abs(s) + s_norm + 1) * rdiff
        return (phi(s+ds, store=False) - phi(s)) / ds

    if search_type == 'wolfe':
        s, phi1, phi0 = nls.scalar_search_wolfe1(phi, derphi, tmp_phi[0],
                                             xtol=1e-2, amin=smin)
    elif search_type == 'armijo':
        s, phi1 = nls.scalar_search_armijo(phi, tmp_phi[0], -tmp_phi[0],
                                       amin=smin)

    if s is None:
        # XXX: No suitable step length found. Take the full Newton step,
        #      and hope for the best.
        s = 1.0

    x = x + s*dx
    if s == tmp_s[0]:
        Fx = tmp_Fx[0]
    else:
        Fx,_,_ = func(x)
    Fx_norm = norm(Fx)

    return s, x, Fx, Fx_norm

def asjacobian(J):
    """
    Convert given object to one suitable for use as a Jacobian.
    """
    if isinstance(J, str):
        return dict(krylov=nls.KrylovJacobian,
                    couplingJacF=CouplingJacobianF,
                    couplingJacCF=CouplingJacobianCF,
                    couplingJacCFT=CouplingJacobianCFT)[J]()   #Laila

class CouplingJacobian(nls.Jacobian):   ## Laila
    r"""
    Find a root of a function, using Krylov approximation for inverse Jacobian.

    This method is suitable for solving large-scale problems.

    Parameters
    ----------
    %(params_basic)s
    rdiff : float, optional
        Relative step size to use in numerical differentiation.
    method : {'lgmres', 'gmres', 'bicgstab', 'cgs', 'minres'} or function
        Krylov method to use to approximate the Jacobian.
        Can be a string, or a function implementing the same interface as
        the iterative solvers in `scipy.sparse.linalg`.

        The default is `scipy.sparse.linalg.lgmres`.
    inner_M : LinearOperator or InverseJacobian
        Preconditioner for the inner Krylov iteration.
        Note that you can use also inverse Jacobians as (adaptive)
        preconditioners. For example,

        >>> from scipy.optimize.nonlin import BroydenFirst, KrylovJacobian
        >>> from scipy.optimize.nonlin import InverseJacobian
        >>> jac = BroydenFirst()
        >>> kjac = KrylovJacobian(inner_M=InverseJacobian(jac))

        If the preconditioner has a method named 'update', it will be called
        as ``update(x, f)`` after each nonlinear step, with ``x`` giving
        the current point, and ``f`` the current function value.
    inner_tol, inner_maxiter, ...
        Parameters to pass on to the \"inner\" Krylov solver.
        See `scipy.sparse.linalg.gmres` for details.
    outer_k : int, optional
        Size of the subspace kept across LGMRES nonlinear iterations.
        See `scipy.sparse.linalg.lgmres` for details.
    %(params_extra)s

    See Also
    --------
    scipy.sparse.linalg.gmres
    scipy.sparse.linalg.lgmres

    Notes
    -----
    This function implements a Newton-Krylov solver. The basic idea is
    to compute the inverse of the Jacobian with an iterative Krylov
    method. These methods require only evaluating the Jacobian-vector
    products, which are conveniently approximated by a finite difference:

    .. math:: J v \approx (f(x + \omega*v/|v|) - f(x)) / \omega

    Due to the use of iterative matrix inverses, these methods can
    deal with large nonlinear problems.

    Scipy's `scipy.sparse.linalg` module offers a selection of Krylov
    solvers to choose from. The default here is `lgmres`, which is a
    variant of restarted GMRES iteration that reuses some of the
    information obtained in the previous Newton steps to invert
    Jacobians in subsequent steps.

    For a review on Newton-Krylov methods, see for example [1]_,
    and for the LGMRES sparse inverse method, see [2]_.

    References
    ----------
    .. [1] D.A. Knoll and D.E. Keyes, J. Comp. Phys. 193, 357 (2004).
           :doi:`10.1016/j.jcp.2003.08.010`
    .. [2] A.H. Baker and E.R. Jessup and T. Manteuffel,
           SIAM J. Matrix Anal. Appl. 26, 962 (2005).
           :doi:`10.1137/S0895479803422014`

    """

    def __init__(self, rdiff=None, method='gmres', inner_maxiter=20,
                 inner_M=None, outer_k=10, **kw):
        self.preconditioner = inner_M
        self.rdiff = rdiff
        self.method = dict(
            bicgstab=scipy.sparse.linalg.bicgstab,
            gmres=scipy.sparse.linalg.gmres,
            lgmres=scipy.sparse.linalg.lgmres,
            cgs=scipy.sparse.linalg.cgs,
            minres=scipy.sparse.linalg.minres,
            ).get(method, method)

        self.method_kw = dict(maxiter=inner_maxiter, M=self.preconditioner)

        if self.method is scipy.sparse.linalg.gmres:
            # Replace GMRES's outer iteration with Newton steps
            self.method_kw['restrt'] = inner_maxiter
            self.method_kw['maxiter'] = 1
            #self.method_kw.setdefault('atol', 0)
        #elif self.method is scipy.sparse.linalg.gcrotmk:
            #self.method_kw.setdefault('atol', 0)
        elif self.method is scipy.sparse.linalg.lgmres:
            self.method_kw['outer_k'] = outer_k
            # Replace LGMRES's outer iteration with Newton steps
            self.method_kw['maxiter'] = 1
            # Carry LGMRES's `outer_v` vectors across nonlinear iterations
            self.method_kw.setdefault('outer_v', [])
            self.method_kw.setdefault('prepend_outer_v', True)
            # But don't carry the corresponding Jacobian*v products, in case
            # the Jacobian changes a lot in the nonlinear step
            #
            # XXX: some trust-region inspired ideas might be more efficient...
            #      See eg. Brown & Saad. But needs to be implemented separately
            #      since it's not an inexact Newton method.
            self.method_kw.setdefault('store_outer_Av', False)
            self.method_kw.setdefault('atol', 0)

        for key, value in kw.items():
            if not key.startswith('inner_'):
                raise ValueError("Unknown parameter %s" % key)
            self.method_kw[key[6:]] = value

    def _update_diff_step(self):
        mx = abs(self.x0).max()
        mf = abs(self.f0).max()
        self.omega = self.rdiff * max(1, mx) / max(1, mf)

    def solve(self, rhs, tol=0):
        if 'tol' in self.method_kw:
            sol, info = self.method(self.op, rhs, **self.method_kw)
        else:
            sol, info = self.method(self.op, rhs, tol=tol, **self.method_kw)
        return sol

    def update(self, x, f):
        self.x0 = x
        self.f0 = f
        self._update_diff_step()

        # Update also the preconditioner, if possible
        if self.preconditioner is not None:
            if hasattr(self.preconditioner, 'update'):
                self.preconditioner.update(x, f)

    def setup(self, x, f, func):
        nls.Jacobian.setup(self, x, f, func)
        self.x0 = x
        self.f0 = f
        self.op = scipy.sparse.linalg.aslinearoperator(self)

        if self.rdiff is None:
            self.rdiff = np.finfo(x.dtype).eps ** (1./2)

        self._update_diff_step()

        # Setup also the preconditioner, if possible
        if self.preconditioner is not None:
            if hasattr(self.preconditioner, 'setup'):
                self.preconditioner.setup(x, f, func)

class CouplingJacobianF(CouplingJacobian):  

    def matvec(self, v):

        nv = norm(v)
        if nv == 0:
            return 0*v
        sc = self.omega / nv
  
        (_, jac_tr_v, jac_ch)= self.func(v)
        
        r = v - np.matmul(jac_ch, v) + np.matmul(jac_ch, jac_tr_v) 
        
        if not np.all(np.isfinite(r)) and np.all(np.isfinite(v)):
            raise ValueError('Function returned non-finite results')
        return r
        
class CouplingJacobianCF(CouplingJacobian):   
  
    def matvec(self, v):

        nv = norm(v)
        if nv == 0:
            return 0*v
        sc = self.omega / nv

        n = v.size
        m=int(n/2)

        dC = v[0:m]
        dF = v[m:n]

        (_, jac_tr_v, jac_ch) = self.func(v)   
	#porosity = 0.97 ### a recuperer des parametres physiques
	#w1 = jac_tr_v*dC+ np.matmul(np.diagflat(porosity), dF)
        w1 = jac_tr_v
        w2 = -np.matmul(jac_ch, dC)+ dF-np.matmul(jac_ch, dF)
        r = np.append(w1, w2)

        if not np.all(np.isfinite(r)) and np.all(np.isfinite(v)):
            raise ValueError('Function returned non-finite results')
        return r
                
class CouplingJacobianCFT(CouplingJacobian):  
  
    def matvec(self, v):

        nv = norm(v)
        if nv == 0:
            return 0*v
        sc = self.omega / nv

        n = v.size
        m=int(n/3)

        dC = v[0:m]
        dT = v[m:2*m]
        dF = v[2*m:n];

        w1 = tr_prime*dC + diag(phi)*dF;
        w2 =  dT - dC - dF;
        w3 = -psi_prime*dT+ dF;

        r = np.append(w1, np.append(w2, w3))  
	
        if not np.all(np.isfinite(r)) and np.all(np.isfinite(v)):
            raise ValueError('Function returned non-finite results')
        return r

  
  
