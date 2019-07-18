import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def func(r,t):
    sigma=10;beta=8/3;rho=28
    x,y,z=r
    return sigma*(y-x),x*(rho-z)-y,x*y - beta*z

def func1(t,r):
    sigma=10;beta=8/3;rho=28
    x,y,z=r.T
    return np.block([[sigma*(y-x)],\
					[x*(rho-z)-(y)],\
					[x*z - beta*z]]).T

n=4000
t=0
dt=0.01
f[0,:] = 1.0
for i in range(n-1):
    f1=dt*func1(t,f[i,:])
    f2=dt*func1(t+(dt/2),f[i,:]+(f1/2))
    f3=dt*func1(t+(dt/2),f[i,:]+(f2/2))
    f4=dt*func1(t+dt,f[i,:]+f3)
    f[i+1,:]=f[i,:]+(f1/6)+(f2/3)+(f3/3)+(f4/6)

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(rs[:,0], rs[:,1], rs[:,2],c='r')
ax.plot(f[:,0], f[:,1], f[:,2],c='b')
plt.show()

t=np.arange(0,40,0.01)
rs=odeint(func,[1,1,1],np.arange(0,40,0.01))

