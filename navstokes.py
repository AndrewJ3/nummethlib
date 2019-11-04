import numpy as np
import scipy.sparse.linalg as spl
import scipy.sparse as sp
from libhelper import relerr
<<<<<<< HEAD
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import sys

=======
#from scipy.fftpack import dct,idct,dst,idst
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import sys
>>>>>>> 43f18a1178fc4c94206979f5514a5e5045a95373
# Timestepping Methods
def rk2(u,v,p,dt,hx,hy,nu,idx,idy,advscheme):
    su1 = advscheme(u,[u,v],idx,idy,hx,hy,dt) + dt*(nu*lap(u,idx,idy,hx,hy) -\
        gradp(p,idx,idy,hx,hy,'x'))
    su2 = (advscheme(u,[u,v],idx,idy,hx,hy,dt) + dt*(nu*lap(u,idx,idy,hx,hy) -\
        gradp(p,idx,idy,hx,hy,'x'))) + dt*su1
    
    sv1 = advscheme(v,[u,v],idx,idy,hx,hy,dt) + dt*(nu*lap(v,idx,idy,hx,hy) -\
            gradp(p,idx,idy,hx,hy,'y'))
    sv2 = (advscheme(v,[u,v],idx,idy,hx,hy,dt) + dt*(nu*lap(v,idx,idy,hx,hy) -\
            gradp(p,idx,idy,hx,hy,'y'))) + dt*sv1
    
    tmpu = u[idx,idy] + (su1 + su2)/2
    tmpv = v[idx,idy] + (sv1 + sv2)/2
    return tmpu,tmpv

def rk4(u,v,p,dt,hx,hy,nu,idx,idy,advscheme):
    su1 = advscheme(u,[u,v],idx,idy,hx,hy,dt) + dt*(nu*lap(u,idx,idy,hx,hy) -\
        gradp(p,idx,idy,hx,hy,'x'))

    su2 = (advscheme(u,[u,v],idx,idy,hx,hy,dt) + dt*(nu*lap(u,idx,idy,hx,hy) -\
        gradp(p,idx,idy,hx,hy,'x'))) + (su1/2)

    su3 = (advscheme(u,[u,v],idx,idy,hx,hy,dt) + dt*(nu*lap(u,idx,idy,hx,hy) -\
        gradp(p,idx,idy,hx,hy,'x'))) + (su2/2)

    su4 = (advscheme(u,[u,v],idx,idy,hx,hy,dt) + dt*(nu*lap(u,idx,idy,hx,hy) -\
        gradp(p,idx,idy,hx,hy,'x'))) + (su3)

    sv1 = advscheme(v,[u,v],idx,idy,hx,hy,dt) + dt*(nu*lap(v,idx,idy,hx,hy) -\
        gradp(p,idx,idy,hx,hy,'y'))

    sv2 = (advscheme(v,[u,v],idx,idy,hx,hy,dt) + dt*(nu*lap(v,idx,idy,hx,hy) -\
        gradp(p,idx,idy,hx,hy,'y'))) + (sv1/2)

    sv3 = (advscheme(v,[u,v],idx,idy,hx,hy,dt) + dt*(nu*lap(v,idx,idy,hx,hy) -\
        gradp(p,idx,idy,hx,hy,'y'))) + (sv2/2)

    sv4 = (advscheme(v,[u,v],idx,idy,hx,hy,dt) + dt*(nu*lap(v,idx,idy,hx,hy) -\
        gradp(p,idx,idy,hx,hy,'y'))) + (sv3)

    tmpu = u[idx,idy]+((su1/6)+(su2/3)+(su3/3)+(su4/6))
    tmpv = v[idx,idy]+((sv1/6)+(sv2/3)+(sv3/3)+(sv4/6))
    return tmpu,tmpv

# advection and diffusion operations
def poifd(f,hy,hx,nx,ny):

    # Sparse Differentiation Matrices (2D Discrete Laplacian)
    l1x=[1]*(nx-1)
    l0x=[-2]*(nx)
    l1y=[1]*(ny-1)
    l0y=[-2]*(ny)
	
    D2_x=sp.diags( [l1x,l0x,l1x], [-1,0,1], format='csr')
    D2_y=sp.diags( [l1y,l0y,l1y], [-1,0,1], format='csr')

    I_x=sp.eye(nx)
    I_y=sp.eye(ny)
    D2x=(hx**-2)*sp.kron(D2_x,I_y)
    D2y=(hy**-2)*sp.kron(I_x,D2_y)
    L=(D2x+D2y)

    #RHS
    f=np.reshape(f,nx*ny)
    
    #krylov methods
<<<<<<< HEAD
    ubicgs,itr=spl.bicgstab(L,f,tol=1e-8)
=======
    ubicgs,itr=spl.bicgstab(L,f,tol=1e-7)
>>>>>>> 43f18a1178fc4c94206979f5514a5e5045a95373
#    print("\r iter = %d"%itr , end = " ")
    
    return ubicgs.reshape(nx,ny)

def lap(u,idx,idy,hx,hy):
	return (u[idx + 1 ,idy] + u[ idx - 1, idy ] - 2*u[idx,idy] )/(hx**2 )+\
	(u[idx , idy + 1 ] + u[ idx , idy - 1 ] -2*u[idx,idy] )/(hy**2)

def upwindord1(u,vel,idx,idy,hx,hy,dt):
    return dt*vel[0][idx,idy]*( u[idx , idy] - u[idx - 1, idy ] )/(hx) +\
               dt*vel[1][idx,idy]*( u[idx , idy  ] - u[idx  , idy - 1])/(hy)

def upwindord2(u,vel,idx,idy,hx,hy,dt):
    return dt*vel[0][idx,idy]*( 3*u[idx , idy] - 4*u[idx - 1, idy ] + u[idx - 2 , idy ])/(2*hx) +\
                   dt*vel[1][idx,idy]*( 3*u[idx , idy] - 4*u[idx , idy - 1 ] + u[idx  , idy - 2 ])/(2*hy)

def laxw(u,vel,idx,idy,hx,hy,dt):
    fin = 0.5*(u[idx + 1, idy] + u[idx,idy]) - (vel[0][idx,idy])*(dt/(2*hx))*(u[idx + 1, idy] - u[idx , idy])
    fout = 0.5*(u[idx , idy] + u[idx - 1,idy]) - (vel[0][idx,idy])*(dt/(2*hx))*(u[idx, idy] - u[idx - 1, idy])
    gin = 0.5*(u[idx + 1, idy] + u[idx,idy]) - (vel[1][idx,idy])*(dt/(2*hy))*(u[idx, idy + 1] - u[idx , idy])
    gout = 0.5*(u[idx , idy] + u[idx,idy - 1]) - (vel[1][idx,idy])*(dt/(2*hy))*(u[idx , idy] - u[idx, idy - 1])
    return (vel[0][idx,idy]*(dt/hx))*( fin - fout ) + (vel[1][idx,idy]*(dt/hy))*( gin - gout )

def ppe(p,u,v,dt,idx,idy,hx,hy,rho):
    nx2 , ny2 = p.shape
    nx = nx2 - 2; ny = ny2 - 2
    
    def gradxu(u,indx,indy,hx):
        return ( u[indx+1,indy] - u[indx-1,indy] )/(2*hx)

    def gradyu(u,indx,indy,hy):
        return ( u[indx,indy+1] - u[indx,indy-1] )/(2*hy)
    
    # rhs \nabla^{2} p = \partial_x u**2  + 2*\partial_y u*\partial_x v + \partial_y u**2
    rhs = np.zeros((nx+2,ny+2))
    rhs[idx,idy] = rho*((hx**2 * hy**2 )/(2*(hx**2 + hy**2 ))) *\
                    ((1/dt)*(gradxu(u,idx,idy,hx) + gradyu(v,idx,idy,hy)) -\
                    ( gradxu(u,idx,idy,hx)**2 + 2*(gradyu(u,idx,idy,hy)*gradxu(v,idx,idy,hx)) +\
                    gradyu(v,idx,idy,hy)**2))

    # solve pressure poisson equation
    p = poifd(rhs[idx,idy],hy,hx,nx,ny)
<<<<<<< HEAD
=======
#     p,res = sorppe(u,v,p,rho,dt,omega,nx,ny,hx,hy)
>>>>>>> 43f18a1178fc4c94206979f5514a5e5045a95373
    return p

#pressure gradient
def gradp(p,idx,idy,hx,hy,component):
	if component == 'x':
		return ( p[ idx + 1 , idy ] - p[ idx - 1 , idy ] )/(2*hx)
	else:
		return ( p[ idx , idy + 1 ] - p[ idx , idy - 1 ] )/(2*hy)

# --------------------------------------------------------------------------
monitorfile = open("../timemonitor","w")
nx = int(sys.argv[1])
ny = int(sys.argv[2])
tfinal = float(sys.argv[3])
dt = float(sys.argv[4])

<<<<<<< HEAD
x0 = -1 ; xn = 3
=======
x0 = -1 ; xn = 5
>>>>>>> 43f18a1178fc4c94206979f5514a5e5045a95373
y0 = -1 ; yn = 1
hy = (yn - y0)/( ny - 1 )
hx = (xn - x0)/( nx - 1 )
u = np.zeros((ny+2,nx+2))
v = np.zeros((ny+2,nx+2))
p = np.zeros((ny+2,nx+2))
xi = np.linspace(x0,xn,nx)
yi = np.linspace(y0,yn,ny)
xx,yy = np.meshgrid(xi,yi)

# constants
nu = 0.1
rho = 1


# interior indices 
uidx = np.arange(2,ny).reshape(ny-2,1)
uidy = np.arange(1,nx+1).reshape(1,nx)
idx = np.arange(1,ny+1).reshape(ny,1)
idy = np.arange(1,nx+1).reshape(1,nx)

# initial conditions 
u[ idx , :1 ] = 1

# BC
u[ -2 ] = 0 
u[ 1  ] = 0 
v[ -2 ] = 0 
v[ 1  ] = 0 
u[ idx , -2 ] = u[ idx , -1  ]
v[ idx , -2 ] = v[ idx , -1  ]
p[ idx , 1  ] = -p[ idx ,  0  ] 
p[ 1 , idy  ] = -p[ 0  , idy  ] 
p[-2 , idy  ] =  p[ -1 , idy  ]
# u[ idx ,-1 ] = (4*u[ idx ,-2 ] - u[ idx ,-3 ])/3
# v[ idx ,-1 ] = (4*v[ idx ,-2 ] - v[ idx ,-3 ])/3
# p[ idx , 0 ] = -(-4*p[ idx , 1 ] + p[ idx , 2 ])/3
# p[ 0 , idy ] = -(-4*p[ 1 , idy ] +  p[ 2 , idy ])/3
# p[-1 , idy ] = (4*p[ -2 , idy ] - p[-3 , idy ])/3
p[ idx , -1 ] =  0
<<<<<<< HEAD
timestep = 0
t = 0
udiff = 1
pdiff = 1
tol = float(sys.argv[5])
while (t < tfinal):
=======
timestep =0
t = 0
udiff = 1
pdiff = 1
while (pdiff > 1e-4 and udiff > 1e-4):
>>>>>>> 43f18a1178fc4c94206979f5514a5e5045a95373
    
    tmpp = ppe(p,u,v,dt,idx,idy,hx,hy,rho)
#     tmpu,tmpv = forwardeuler(u,v,dt,hx,hy,nu,idx,idy,laxw)
    tmpu,tmpv = rk4(u,v,p,dt,hx,hy,nu,idx,idy,laxw)
#     tmpu,tmpv = rk2(u,v,dt,hx,hy,nu,idx,idy,upwindord1)
    
    # update boundary conditions
    tmpu[ : , -1 ] = tmpu[ : , -2  ]
    tmpv[ : , -1 ] = tmpv[ : , -2  ]
    tmpp[ : , 0  ] = tmpp[ : , 1  ] 
    tmpp[ 0 , :  ] = tmpp[ 1 , :  ] 
    tmpp[ -1 , :  ] = tmpp[ -2 , :  ]   
    tmpp[ : , -1 ] = 0
    tmpu[ -1 ] = 0 
    tmpu[ 0 ] = 0 
    tmpv[ -1 ] = 0 
    tmpv[ 0 ] = 0 
    
    # Momentum and Pressure Tolerance
    udiff = relerr(tmpu,u[idx,idy]+2.3e-16)
    pdiff = relerr(tmpp,p[idx,idy]+2.3e-16)
    # update momentum and pressure for next timestep
    u[idx,idy] = tmpu
    v[idx,idy] = tmpv
    p[idx,idy] = tmpp
	
    t += dt
    #print("\r time = %g"%t , end = " ")
	# Advance and Display/Monitor Time Step
    t = round(t,3)
<<<<<<< HEAD
=======
    pdiff = round(pdiff,4)
    udiff = round(udiff,4)
>>>>>>> 43f18a1178fc4c94206979f5514a5e5045a95373
    timestep += 1
    monitorfile.write(" Time: " + str(t) + " Number of Timesteps: "+ str(timestep)+"\n")
    monitorfile.write("momentum tolerance: " + str(udiff)+"\n")
    monitorfile.write("pressure tolerance: "+ str(pdiff)+"\n")
<<<<<<< HEAD
    monitorfile.write("------------------------------------\n")
    if pdiff < tol and udiff < tol:
          break;
=======
    monitorfile.write("------------------------\n")
>>>>>>> 43f18a1178fc4c94206979f5514a5e5045a95373
#print(" ")

#zip data
np.savez('test.npz',u=u,v=v,p=p,x=xx,y=yy)
<<<<<<< HEAD
um = np.sqrt(u**2 + v**2)
print(um[idx,-2].max())
ums = um[idx,-2]/np.max(um[idx,-2],0)
ums = ums.flatten()

#visualization
fig = plt.figure(figsize=(14,8))
ax = fig.add_subplot(1,1,1)
plt1 = ax.contourf(xx,yy,u[idx,idy],100,cmap='jet')
=======

#visualization
xx = xx[idx,idy]; yy = yy[idx,idy]
fig = plt.figure(figsize=(14,8))
ax = fig.add_subplot(1,1,1)
plt1 = ax.contourf(xx,yy,np.sqrt(u[idx,idy]**2 + v[idx,idy]**2),100,cmap='jet')
>>>>>>> 43f18a1178fc4c94206979f5514a5e5045a95373
fig.colorbar(plt1)
plt.savefig('umag.png')

fig = plt.figure(figsize=(14,8))
ax = fig.add_subplot(1,1,1)
plt2 = ax.contourf(xx,yy,p[idx,idy],100,cmap='jet')
fig.colorbar(plt2)
plt.savefig('press.png')
<<<<<<< HEAD

fig = plt.figure(figsize=(14,8))
ax = fig.add_subplot(1,1,1)
ax.plot(yi,1-yi**2,'b--',label = 'exact')
ax.plot(yi,ums,'ro',label = 'FD')
plt.legend()
plt.savefig('error.png')

=======
>>>>>>> 43f18a1178fc4c94206979f5514a5e5045a95373
