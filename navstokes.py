def lap(u,idx,idy,hx,hy):
	return (u[idx + 1 ,idy] + u[ idx - 1, idy ] - 2*u[idx,idy] )/(hx**2 )+\
	(u[idx , idy + 1 ] + u[ idx , idy - 1 ] -2*u[idx,idy] )/(hy**2)

def upwind(u,vel,idx,idy,hx,hy):
	return vel[0][idx,idy]*( u[idx , idy] - u[idx - 1, idy ] )/hx +\
               vel[1][idx,idy]*( u[idx , idy] - u[idx , idy - 1] )/hy

# def gradxu(u,indx,indy,hx):
#     return ( u[indx+1,indy] - u[indx-1,indy] )/(2*hx)

# def gradyu(u,indx,indy,hy):
#     return ( u[indx,indy+1] - u[indx,indy-1] )/(2*hy)

def ppe(p,u,v,dt,idx,idy,hx,hy,rho):
    nx2 , ny2 = p.shape
    nx = nx2 - 2; ny = ny2 - 2
    
    def gradxu(u,indx,indy,hx):
        return ( u[indx+1,indy] - u[indx-1,indy] )/(2*hx)

    def gradyu(u,indx,indy,hy):
        return ( u[indx,indy+1] - u[indx,indy-1] )/(2*hy)
    
    omega = 0.99
    tol = 1e-3
    maxiter = 1000
    for ii in range(maxiter):
        
        p[idx,idy] = (0.25*omega)*(((p[idx+1,idy] + p[idx-1,idy])*hy**2 + (p[idx,idy+1] + p[idx,idy-1])*hx**2)/(2*(hx**2 + hy**2)) -\
            rho*((hx**2 * hy**2 )/(2*(hx**2 + hy**2 ))) *((1/dt)*(gradxu(u,idx,idy,hx) + gradyu(v,idx,idy,hy)) -\
            ( gradxu(u,idx,idy,hx)**2 + 2*(gradyu(u,idx,idy,hy)*gradxu(v,idx,idy,hx)) + gradyu(v,idx,idy,hy)**2))) - (1.0-omega)*p[idx,idy]


        res = np.zeros((nx+1,ny+1))
        res[idx,idy] = np.abs(p[idx,idy] - (((p[idx+1,idy] + p[idx-1,idy])*hy**2 + (p[idx,idy+1] + p[idx,idy-1])*hx**2)/(2*(hx**2 + hy**2)) -\
            rho*((hx**2 * hy**2 )/(2*(hx**2 + hy**2 ))) *((1/dt)*(gradxu(u,idx,idy,hx) + gradyu(v,idx,idy,hy)) -\
            ( gradxu(u,idx,idy,hx)**2 + 2*(gradyu(u,idx,idy,hy)*gradxu(v,idx,idy,hx)) + gradyu(v,idx,idy,hy)**2))))

        
#         for i in range(1,nx+1):
#             for j in range(1,ny+1):
#                 p[i,j] = ((p[i+1,j] + p[i-1,j])*hy**2 + (p[i,j+1] + p[i,j-1])*hx**2)/(2*(hx**2 + hy**2)) -\
#                 rho*((hx**2 * hy**2 )/(2*(hx**2 + hy**2 ))) *((1/dt)*(gradxu(u,i,j,hx) + gradyu(v,i,j,hy)) -\
#                 ( gradxu(u,hx,i,j)**2 + 2*(gradyu(u,hy,i,j)*gradxu(v,hx,i,j)) + gradyu(v,hy,i,j)**2))


#         res = np.zeros((nx+1,ny+1))
#         for i in range(1,nx):
#             for j in range(1,ny):
#                 res[i,j] = np.abs(p[i,j] - ((p[i+1,j] + p[i-1,j])*hy**2 + (p[i,j+1] + p[i,j-1])*hx**2)/(2*(hx**2 + hy**2)) -\
#                 rho*((hx**2 * hy**2 )/(2*(hx**2 + hy**2 ))) *((1/dt)*(gradxu(u,i,j,hx) + gradyu(v,i,j,hy)) -\
#                 ( gradxu(u,hx,i,j)**2 + 2*(gradyu(u,hy,i,j)*gradxu(v,hx,i,j)) + gradyu(v,hy,i,j)**2)))
        
        if tol <= norm(res):
            break;

    return p


def gradp(p,idx,idy,hx,hy,component):
	if component == 'x':
		return ( p[ idx + 1 , idy ] - p[ idx - 1 , idy ] )/(2*hx)
	else:
		return ( p[ idx , idy + 1 ] - p[ idx , idy - 1 ] )/(2*hy)

#def laxfr(u,vel,idx,idy,h):
#	return 0.5*(u[idx+1,idy] - u[idx-1,idy] ) +\
#               0.5*vel[0][idx,idy]*( u[idx + 1 , idy] - u[idx - 1, idy ] )/h +\
#               0.5*vel[1][idx,idy]*( u[idx , idy + 1] - u[idx , idy - 1] )/h
nx = 50//2
ny = 100//2
tfinal = 10
dt = 0.01

x0 = -1 ; xn = 2
y0 = -1 ; yn = 1
hy = (yn - y0)/( nx + 1 )
hx = (xn - x0)/( ny + 1 )
u = np.zeros((nx+2,ny+2))
v = np.zeros((nx+2,ny+2))
p = np.zeros((nx+2,ny+2))
yi = np.linspace(y0,yn,nx+2)
xi = np.linspace(x0,xn,ny+2)
xx,yy = np.meshgrid(xi,yi)

# constants
nu = 0.1
rho = 1

# interior indices 
idx = np.arange(1,nx+1).reshape(nx,1)
idy = np.arange(1,ny+1).reshape(1,ny)

# initial/ boundary conditions
u[ : , 0 ] = 1
u[ -1 , :] = 0 
u[ 0 , :] = 0 
u[ : ,-1] = -(-4*u[ : , -2] + u[ : ,-3])/3
p[ : , -1 ] = -(-4*p[ : , -2 ] + p[ : , -3 ])/3
p[ : , -1 ] = 0
p[ 0 , : ] = -(-4*p[ 1 , : ] +  p[ 2 , : ])/3
p[-1 , : ] = -(-4*p[ -2 , : ] + p[-3 , : ])/3

t = 0
while (t < tfinal):

    # solve momentum and pressure equations
    tmpu = u[idx,idy] + dt*(-upwind(u,[u,v],idx,idy,hx,hy)+nu*lap(u,idx,idy,hx,hy) -\
            gradp(p,idx,idy,hx,hy,'x'))
    tmpv = v[idx,idy] + dt*(-upwind(v,[u,v],idx,idy,hx,hy)+nu*lap(v,idx,idy,hx,hy) -\
            gradp(p,idx,idy,hx,hy,'y'))
    tmpp = ppe(p,u,v,dt,idx,idy,hx,hy,rho)[idx,idy]

    # update boundary conditions
    u[ : ,-1 ] = -(-4*u[ : ,-2] + u[ : ,-3])/3
    u[ -1 , : ] = 0 
    u[ 0 , : ] = 0 
    p[ : , 0 ] = -(-4*p[ : , -2 ] + p[ : , -3 ])/3
    p[ : , -1 ] = 0
    p[ 0 , : ] = -(-4*p[ 1 , : ] + p[ 2 , : ])/3
    p[-1 , : ] = -(-4*p[ -2 , : ] + p[-3 , : ])/3


    # update momentum and pressure for next timestep
    u[idx,idy] = tmpu
    v[idx,idy] = tmpv
    p[idx,idy] = tmpp

    t += dt
    print("\r time = %g"%t , end = " ")

print(" ")
#print(u)
#visualization
fig=plt.figure(figsize=(20,12))
ax=fig.add_subplot(1,1,1)
plt1 = ax.contourf(xx,yy,np.sqrt(u**2+v**2),100,cmap ='viridis')
fig.colorbar(plt1)
# plt.savefig('nse.png')

