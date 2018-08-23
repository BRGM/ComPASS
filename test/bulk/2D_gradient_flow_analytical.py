import numpy as np
import matplotlib.pyplot as plt
import time
import math
from scipy.special import erfc


#calcul de la solution exacte
Phi = 0.2
rhow = 1000

b = 5e+03
rhoCp = 1000*2000
Cm = 60+273   #température pour x=0
gamma = 0
mu = 0
U = 1/(Phi+rhoCp*(1-Phi)/(rhow*b))*b #*densité molaire à priori constante
Dx = 1
Rd = rhow*Phi*b+rhoCp*(1-Phi)
P = (U**2 + 4*Dx*Rd*(mu+gamma))**(1/2)

def exact_sol(x,t):
    x = np.asarray(x)
    _where = x < 200
    exact_sol = np.zeros(x.shape, dtype=np.double)
    xw = x[_where]
    type(xw)
    print(Dx,Rd,t)
    erfc((Rd*xw-P*t)/(2*(Dx*Rd*t)**0.5))
    type(exact_sol[_where])
    exact_sol[_where,t] = 1/2*Cm*np.exp(gamma*t)*(np.exp(xw*(U-P)/(2*Dx))*erfc((Rd*xw-P*t)/(2*(Dx*Rd*t)**0.5))+np.exp(xw*(U-P)/(2*Dx))*erfc((Rd*xw-P*t)/(2*(Dx*Rd*t)**0.5)))
    #exact_sol[_where] = 1.
    #where = np.logical_not(_where)
    #xw = x[_where]
    return exact_sol

Lx, Ly, Lz = 100., 10., 1.
dx = 2

day = 2. # seconds choix!!!!! relation avec np.arrange()
x = np.arange(0, Lx , dx)
t = day * np.hstack([0.1, np.arange(20, 1020, 30), 1E8/day])
table = np.array([exact_sol(ti, x) for ti in t>0])

csv_table = np.hstack([np.reshape(t/day, (-1, 1)),np.reshape(t, (-1, 1)),table])

np.savetxt('2D_GF_ANALYTICAL.csv', csv_table,
        delimiter=';',
        header='seconds;' + ';'.join('%.1f'%xi for xi in x)
)

try:
    import matplotlib.pyplot as plt
except ImportError:
    print('No figure is drawn as matpplotlib is not installed.')
    raise
else:
    plt.clf()
    plt.plot(x, table[-1])
    plt.xlabel('distance (m)')
    plt.ylabel('temperature (deg C)')
    plt.grid(True)
    plt.savefig('2D_GF_ANALYTICAL.png')
