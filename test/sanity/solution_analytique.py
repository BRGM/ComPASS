import numpy as np
import matplotlib.pyplot as plt
from scipy import special
import time
import math

"""def erfterm(x, n):

def erf(x):
  def term(n):
      num = (-1.0)**n * x**(2*n+1)
      denum = fact(n) * (2*n+1)
      return num/denum
   l = sum((erfterm(x, i) for i in range(25)))*(2*sqrt(np.pi))
   return l"""

def f(tau,x,y):
    return tau**(-3/2)*(special.erf((a+y)/math.sqrt(4*Dt*tau))+special.erf((a-y)/math.sqrt(4*Dt*tau)))*np.exp(-(x-tau*v)**2/(4*Dl*tau)) #y pas encore défini


def exact_sol(x,y,t):
    n=100
    s=0
    for i in range(n):
        tau=t*np.cos((2*(i+1)-1)*np.pi/(4*n))**2
        s=s+x*t*c0/(32*np.pi*Dl)**2 * np.pi/n*np.sin((2*(i+1)-1)*np.pi/(2*n))*f(tau,x,y)
    return s


"""for i in range(nx):
    tab[i] = ny*[0]
for i in range(nx):
    for j in range(ny):
        tab[i][j] = 2*[0]
for i in range(nx):
    for j in range(ny):
        tab[i][j][0]= dx*i
        tab[i][j][1]= dy*j"""


def test():
    nx = math.floor(Lx/dx)
    ny = math.floor(Ly/dy)
    x=np.linspace(0,Lx,nx)
    y=np.linspace(-Ly/2,Ly/2,ny)
    xx, yy = np.meshgrid(x, y)
    
    for k in range(nt):
        t=(k+1)*dt
        c=exact_sol(xx,yy,t)
        if ((k+1)%10==0):
            h = plt.contourf(x,y,c)
            plt.title('t='+str(t))
            plt.show()
            #time.sleep(10)

if __name__ == '__main__':
    Dl = 1
    Dt = 0.1
    c0 = 10
    a = 8
    Pe = 10
    v = 1
    Lx , Ly = 100 , 40
    dx , dy = 2 , 2  #dx=xy=2
    tfinal = 100
    dt = 0.2
    nt = math.floor(tfinal/dt)
    test()
