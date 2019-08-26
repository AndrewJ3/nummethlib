import numpy as np
import scipy.sparse.linalg as spl
import scipy.sparse as sp
#from scipy.fftpack import dct,idct,dst,idst
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import sys
# Timestepping Methods
def rk2(u,v,dt,hx,hy,nu,idx,idy,advscheme):
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

def rk4(u,v,dt,hx,hy,nu,idx,idy,advscheme):
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
def poibicgstab(f,hy,hx,nx,ny):

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
    ubicgs,itr=spl.bicgstab(L,f,tol=1e-7)
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
    p = poibicgstab(rhs[idx,idy],hy,hx,nx,ny)
#     p,res = sorppe(u,v,p,rho,dt,omega,nx,ny,hx,hy)
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

x0 = -1 ; xn = 3
y0 = -1 ; yn = 1
hy = (yn - y0)/( ny + 2 )
hx = (xn - x0)/( nx + 2 )
u = np.zeros((ny+2,nx+2))
v = np.zeros((ny+2,nx+2))
p = np.zeros((ny+2,nx+2))
xi = np.linspace(x0,xn,nx+2)
yi = np.linspace(y0,yn,ny+2)
xx,yy = np.meshgrid(xi,yi)

# constants
nu = 0.1
rho = 1

# interior indices 
idx = np.arange(1,ny+1).reshape(ny,1)
idy = np.arange(1,nx+1).reshape(1,nx)

# initial conditions 
u[ : , 0 ] = 1

# boundary conditions
u[ -1 , : ] = 0 
u[ 0 , : ] = 0 
v[ -1 , : ] = 0 
v[ 0 , : ] = 0 
# u[ : , -1 ] = u[ : , -2 ]
# v[ : , -1 ] = v[ : , -2 ]
u[ : ,-1 ] = (4*u[ : ,-2 ] - u[ : ,-3 ])/3
v[ : ,-1 ] = (4*v[ : ,-2 ] - v[ : ,-3 ])/3
p[ : , 0 ] = -(-4*p[ : , 1 ] + p[ : , 2 ])/3
p[ 0 , : ] = -(-4*p[ 1 , : ] +  p[ 2 , : ])/3
p[-1 , : ] = (4*p[ -2 , : ] - p[-3 , : ])/3
# p[ : , 0 ] = -p[ : , 1 ]
# p[ 0 , : ] = -p[ 1 , : ] 
# p[-1 , : ] =  p[ -2 , : ] 
p[ : , -1 ] = 0

t = 0
timestep = 0
while (t < tfinal):
    
    tmpp = ppe(p,u,v,dt,idx,idy,hx,hy,rho)
    tmpu,tmpv = rk4(u,v,dt,hx,hy,nu,idx,idy,laxw)
    
    # update boundary conditions
    u[ -1 , : ] = 0 
    u[ 0 , : ] = 0 
    v[ -1 , : ] = 0 
    v[ 0 , : ] = 0 
#     v[ : , -1 ] = v[ : , -2 ]
#     u[ : , -1 ] = u[ : , -2 ]
    u[ : ,-1 ] = (4*u[ : ,-2 ] - u[ : ,-3 ])/3
    v[ : ,-1 ] = (4*v[ : ,-2 ] - v[ : ,-3 ])/3
#     p[ : , 0 ] = -p[ : , 1 ]
#     p[ 0 , : ] = -p[ 1 , : ] 
#     p[-1 , : ] =  p[ -2 , : ]
    p[ : , 0 ] = -(-4*p[ : , 1 ] + p[ : , 2 ])/3
    p[ 0 , : ] = -(-4*p[ 1 , : ] +  p[ 2 , : ])/3
    p[-1 , : ] = (4*p[ -2 , : ] - p[-3 , : ])/3
    p[ : , -1 ] = 0
    
#     [idx,idy]

    # update momentum and pressure for next timestep
    u[idx,idy] = tmpu
    v[idx,idy] = tmpv
    p[idx,idy] = tmpp

    t += dt
    #print("\r time = %g"%t , end = " ")
	# Advance and Display/Monitor Time Step
    t = round(t,3)
    timestep += 1
    monitorfile.write(" Time: " + str(t) + " Number of Timesteps: "+ str(timestep)+"\n")
    monitorfile.write("------------------------\n")
#print(" ")

#zip data
np.savez('test.npz',u= u,v=v,p=p,x=xx,y=yy)

#visualization
fig = plt.figure(figsize=(14,8))
ax = fig.add_subplot(1,1,1)
plt1 = ax.contourf(xx,yy,np.sqrt(u**2 + v**2),100,cmap='jet')
fig.colorbar(plt1)
plt.savefig('umag.png')

fig = plt.figure(figsize=(14,8))
ax = fig.add_subplot(1,1,1)
plt2 = ax.contourf(xx,yy,p,100,cmap='jet')
fig.colorbar(plt2)
plt.savefig('press.png')
