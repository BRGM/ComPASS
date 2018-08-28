import numpy as np
import matplotlib.pyplot as plt
import time
import math
from math import *
from scipy.special import erfc


Phi = 0.2
rhow = 1000
b = 5e+03
rhoCp = 1000*2000
Rd = rhow*Phi*b+rhoCp*(1-Phi)
mu = b*rhow/Rd
C0 = 60+273   #température pour x=0
gamma = 0 #pas sur
mu = 1
Lx, Ly, Lz = 100., 50., 1.
n0 = Ly
epsilon = 8
n = 100
#U = 1/(Phi+rhoCp*(1-Phi)/(rhow*b))*b #*densité molaire à priori constante
#Dx = 1
#P = (U**2 + 4*Dx*Rd*(mu+gamma))**(1/2)

def exact_sol(x,y,t):
    #exact_sol = np.zeros(x.shape, dtype=np.double)
    print(Dx,Rd,t)
    exact_sol(x,y,t) = epsilon*C0/(2*n0)*(np.erfc((x-mu*t)/sqrt(4*Dl*t))+np.exp(mu*x/Dl)*np.erfc((x+mu*t)/sqrt(4*Dl*t))
    for i in range(n)
        l1 = C0/(i*np.pi)*np.sin(i*np.pi*epsilon/n0)*np.cos(i*np.pi*yw/n0)
        l2 = np.exp(x/2*(mu/Dl-sqrt(mu/Dl+(2*i*np.pi/n0)**2*Dt/Dl))
        l3 = np.erfc((x-Dl*t*((mu/Dl)**2+(2*i*np.pi/n0)**2*Dt/Dl)/sqrt(4*Dl*t)))
        l4 = np.exp(x/2*(mu/Dl+sqrt((mu/Dl)**2+(2*i*np.pi/n0)**2*Dt/Dl)))
        l5 = np.erfc((x+Dl*t*sqrt((mu/Dl)**2+(2*i*np.pi/n0)*Dt/Dl))/sqrt(4*Dl*t))
        exact_sol(x,y,t) = exact_sol(x,y,t) + l1-l2-l3+l4-l5
    return exact_sol

t = 100
dt = 2
nt = math.floor(t/dt)

def test():
    nx = 100
    ny = 50
    x=np.linspace(0,Lx,nx)
    y=np.linspace(-Ly/2,Ly/2,ny)
    xx, yy = np.meshgrid(x, y)

    for k in range(nt):
        c=exact_sol(xx,yy,t)
        if ((k+1)%10==0):
            plt.subplot(211)
            cs = plt.contourf(x,y,c2)
            fig.colorbar(cs)
            plt.draw()
            plt.pause(0.1)
            plt.clf()
