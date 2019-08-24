import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt
import sys
from mpl_toolkits import mplot3d

def lap(u,idx,idy,hx,hy):
	return (u[idx + 1 ,idy] + u[idx + 1,idy] - 2*u[idx,idy])/hx**2 +\
	(u[idx , idy + 1 ] + u[idx , idy - 1] -2*u[idx,idy])/hy**2

def upwind(u,vel,idx,idy,hx,hy):
	return vel[0][idx,idy]*( u[idx , idy] - u[idx - 1, idy ] )/hx +\
               vel[1][idx,idy]*( u[idx , idy] - u[idx , idy - 1] )/hy

def ppe(p,u,v,dt,idx,idy,h,rho):

	def gradxu(u,h,i,j):
		return ( u[i+1,j] - u[i-1,j] )/(2*h)

	def gradyu(u,h,i,j):
		return ( u[i,j+1] - u[i,j-1] )/(2*h)
	
	tol = 1e-4
	for i in range(1,n-1):
		for j in range(1,n-1):
			p[i,j] = ( p[i+1,j] + p[i-1,j] + p[i,j-1] + p[i,j+1] ) +\
            4*(h**2)*rho*( gradxu(u,


	res = np.zeros((n-1,n-1))
	for i in range(1,n-1):
		for j in range(1,n-1):
			res[ i , j ] = p[i,j] - (h**2) * ( p[idx + 1 , idy] + p[idx - 1, idy ] +\
								 p[idx - 1, idy ] + p[ idx - 1, idy ] )

	
	return p[idx,idy]

def pressgrad(p,idx,idy,hx,hy,component):
	if component == 'x':
		return ( p[ idx + 1 , idy ] - p[ idx - 1 , idy ] )/(2*h)
	else:
		return ( p[ idx , idy + 1 ] - p[ idx , idy - 1 ] )/(2*h)

#def laxfr(u,vel,idx,idy,h):
#	return 0.5*(u[idx+1,idy] - u[idx-1,idy] ) +\
#               0.5*vel[0][idx,idy]*( u[idx + 1 , idy] - u[idx - 1, idy ] )/h +\
#               0.5*vel[1][idx,idy]*( u[idx , idy + 1] - u[idx , idy - 1] )/h

n = int(sys.argv[1])
tfinal = float(sys.argv[2])
dt = float(sys.argv[3])

x0 = -1 ; xn = 2
y0 = -1 ; yn = 1
hy = (yn - y0)/( n + 1 )
hx = (xn - x0)/( n + 1 )
u = np.zeros((n+2,n+2))
v = np.zeros((n+2,n+2))
p = np.zeros((n+2,n+2))
yi = np.linspace(y0,yn,n)
xi = np.linspace(x0,xn,n)
xx,yy = np.meshgrid(xi,yi)

# constants
nu = 0.1
rho = 1

# initial conditions
u[ : , 0 ] = 1

# interior indices 
idx = np.arange(1,n+1).reshape(n,1)
idy = idx.T

t = 0
while (t <= tfinal):
	tmpu = u[idx,idy] + dt*(-upwind(u,[u,v,],idx,idy,h)+nu*lap(u,idx,idy,h) -\
			pressgrad(p,idx,idy,h,'x'))
	tmpv = v[idx,idy] + dt*(-upwind(v,[u,v],idx,idy,h)+nu*lap(v,idx,idy,h) -\
			pressgrad(p,idx,idy,h,'y'))
	u[idx,idy] = tmpu
	v[idx,idy] = tmpv
	t += dt
	print("\r time = %g"%t , end = " ")
	
print(" ")

#visualization
fig=plt.figure(figsize=(20,12))
ax=fig.add_subplot(1,1,1)
plt1 = ax.contourf(xx,yy,u[idx,idy],100,cmap ='viridis')
fig.colorbar(plt1)
plt.savefig('nse.png')


