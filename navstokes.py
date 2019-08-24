def lap(u,idx,idy,hx,hy):
	return (u[idx + 1 ,idy] + u[ idx + 1, idy ] - 2*u[idx,idy] )/hx**2 +\
	(u[idx , idy + 1 ] + u[ idx , idy - 1 ] -2*u[idx,idy] )/hy**2

def upwind(u,vel,idx,idy,hx,hy):
	return vel[0][idx,idy]*( u[idx , idy] - u[idx - 1, idy ] )/hx +\
               vel[1][idx,idy]*( u[idx , idy] - u[idx , idy - 1] )/hy

#def lap(u,hx,hy):
#	return (u[: - 1 ,:] + u[ : + 1, : ] - 2*u[:,:] )/hx**2 +\
#			(u[ : , : + 1 ] + u[ : , : - 1 ] -2*u[:,:] )/hy**2

#def upwind(u,vel,hx,hy):
#	return vel[0]*( u[ : , : ] - u[ : - 1, : ] )/hx +\
#               vel[1]*( u[ : , : ] - u[ : , : - 1] )/hy

def ppe(p,u,v,dt,idx,idy,hx,hy,rho):

	def gradxu(u,hx,i,j):
		return ( u[i+1,j] - u[i-1,j] )/(2*hx)

	def gradyu(u,hy,i,j):
		return ( u[i,j+1] - u[i,j-1] )/(2*hy)
	
	tol = 1e-3
    maxiter = 1000
    p[idx,idy]
    for ii in range(maxiter)
        for i in range(1,n+1):
            for j in range(1,n+1):
                p[i,j] = ((p[i+1,j] + p[i-1,j])*hy**2 + (p[i,j+1] + p[i,j-1])*hx**2)/(2*(hx**2 + hy**2)) -\
                rho*((hx**2 * hy**2 )/(2*(hx**2 + hy**2 ))) *((1/dt)*(gradxu(u,i,j,hx) + gradyu(v,i,j,hy)) -\
               ( gradxu(u,hx,i,j)**2 + 2*(gradyu(u,hy,i,j)*gradxu(v,hx,i,j)) + gradyu(v,hy,i,j)**2))


        res = np.zeros((n+1,n+1))
        for i in range(1,n):
            for j in range(1,n):
                res[i,j] = np.abs(p[i,j] - ((p[i+1,j] + p[i-1,j])*hy**2 + (p[i,j+1] + p[i,j-1])*hx**2)/(2*(hx**2 + hy**2)) -\
                rho*((hx**2 * hy**2 )/(2*(hx**2 + hy**2 ))) *((1/dt)*(gradxu(u,i,j,hx) + gradyu(v,i,j,hy)) -\
               ( gradxu(u,hx,i,j)**2 + 2*(gradyu(u,hy,i,j)*gradxu(v,hx,i,j)) + gradyu(v,hy,i,j)**2)))
        
        if tol <= norm(res):
            break;

	return p


def pressgrad(p,idx,idy,hx,hy,component):
	if component == 'x':
		return ( p[ idx + 1 , idy ] - p[ idx - 1 , idy ] )/(2*hx)
	else:
		return ( p[ idx , idy + 1 ] - p[ idx , idy - 1 ] )/(2*hy)

#def pressgrad(p,hx,hy,component):
#	if component == 'x':
#		return ( p[ : + 1 , : ] - p[ : - 1 , : ] )/(2*hx)
#	else:
#		return ( p[ : , : + 1 ] - p[ : , : - 1 ] )/(2*hy)

#def laxfr(u,vel,idx,idy,h):
#	return 0.5*(u[idx+1,idy] - u[idx-1,idy] ) +\
#               0.5*vel[0][idx,idy]*( u[idx + 1 , idy] - u[idx - 1, idy ] )/h +\
#               0.5*vel[1][idx,idy]*( u[idx , idy + 1] - u[idx , idy - 1] )/h

n = 25
tfinal = 100
dt = 0.01

x0 = -1 ; xn = 2
y0 = -1 ; yn = 1
hy = (yn - y0)/( n + 1 )
hx = (xn - x0)/( n + 1 )
u = np.zeros((n+2,n+2))
v = np.zeros((n+2,n+2))
p = np.zeros((n+2,n+2))
yi = np.linspace(y0,yn,n+2)
xi = np.linspace(x0,xn,n+2)
xx,yy = np.meshgrid(xi,yi)

# constants
nu = 0.1
rho = 1

# interior indices 
idx = np.arange(1,n+1).reshape(n,1)
idy = idx.T
print(len(idx))
# initial conditions
u[ idx , 0 ] = 1

t = 0
while (t < tfinal):
	tmpu = u[idx,idy] + dt*(-upwind(u,[u,v],idx,idy,hx,hy)+nu*lap(u,idx,idy,hx,hy) -\
			pressgrad(p,idx,idy,hx,hy,'x'))
	tmpv = v[idx,idy] + dt*(-upwind(v,[u,v],idx,idy,hx,hy)+nu*lap(v,idx,idy,hx,hy) -\
			pressgrad(p,idx,idy,hx,hy,'y'))
	u[idx,-1] = (-4*u[idx ,-2] + u[idx , -3 ])/3
	u[ -1 , :] = 0 
	u[ 0 , :] = 0 

	u[idx,idy] = tmpu
	v[idx,idy] = tmpv

#while (t < tfinal):
#	tmpu = u + dt*(-upwind(u,[u,v],hx,hy)+nu*lap(u,hx,hy) -\
#			pressgrad(p,hx,hy,'x'))
#	tmpv = v + dt*(-upwind(v,[u,v],hx,hy)+nu*lap(v,hx,hy) -\
#			pressgrad(p,hx,hy,'y'))
#	u[:,-1] = (-4*u[: ,-2] + u[: , -3 ])/3
#	u[ -1 , :] = 0 
#	u[ 0 , :] = 0 
#	u = tmpu
#	v = tmpv


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

