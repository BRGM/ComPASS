import numpy as np
import matplotlib.pyplot as plt
import ComPASS
import time
import math
from math import *
from scipy import special



#U = 1/(Phi+rhoCp*(1-Phi)/(rhow*b))*b #*densité molaire à priori constante
#Dx = 1
#P = (U**2 + 4*Dx*Rd*(mu+gamma))**(1/2)

def exact_sol(x,y,t):
    exact_sol = epsilon*C0/(2*n0)*(special.erfc((x-mu*t)/sqrt(4*Dl*t))+np.exp(mu*x/Dl)*special.erfc((x+mu*t)/sqrt(4*Dl*t)))
    for i in range(1,n):
        l1 = C0/(i*np.pi)*np.sin(i*np.pi*epsilon/n0)*np.cos(i*np.pi*y/n0)
        l2 = np.exp(x/2*(mu/Dl-sqrt((mu/Dl)**2+(2*i*np.pi/n0)**2*Dt/Dl)))
        l3 = special.erfc((x-Dl*t*sqrt((mu/Dl)**2+(2*i*np.pi/n0)**2*Dt/Dl))/sqrt(4*Dl*t))
        l4 = np.exp(x/2*(mu/Dl+sqrt((mu/Dl)**2+(2*i*np.pi/n0)**2*Dt/Dl)))
        l5 = special.erfc((x+Dl*t*sqrt((mu/Dl)**2+(2*i*np.pi/n0)**2*Dt/Dl))/sqrt(4*Dl*t))
        exact_sol = exact_sol +l1*(l2*l3+l4*l5)
    return exact_sol



def test():
    nx = 100
    ny = 50
    nz = 1
    x=np.linspace(0,Lx,nx)
    y=np.linspace(-Ly/2,Ly/2,ny)
    xx, yy = np.meshgrid(x, y)
    for k in range(nt):
        t = (k+1)*dt
        print(k)
        c=exact_sol(xx,yy,t)
        if ((k+1)%5==0):
            fig = plt.figure(1)
            plt.subplot(211)
            cs = plt.contourf(x,y,c)
            fig.colorbar(cs)
            plt.title('t='+str(t))
            plt.draw()
            plt.pause(0.1)
            #plt.savefig(('exact_sol_'+str(t)),format='png')
            plt.clf()


if __name__ == '__main__':
    tfinal = 100
    dt = 2
    nt = math.floor(tfinal/dt)
    Phi = 0.2
    rhow = 1000
    b = 5e+03
    rhoCp = 1000*2000
    Rd = rhow*Phi*b+rhoCp*(1-Phi)
    print(Rd)
    mu = b*rhow/Rd
    C0 = 60+273  #température en x=0
    gamma = 0 #pas sur dans l'expression de P
    mu = 1
    Dl = 1e+7/Rd  #conductivité thermique/Rd
    Dt = 1e+7/Rd
    print(Dt)
    Lx, Ly, Lz = 100., 50., 1.
    n0 = Ly
    epsilon = 8
    n = 100
    test()
