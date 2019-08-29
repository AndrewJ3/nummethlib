import numpy as np
from libhelper import relerr
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

def lap(u,idx,idy,hx,hy):
    return (u[idx + 1 ,idy] + u[ idx - 1, idy ] - 2*u[idx,idy] )/(hx**2 )+\
    (u[idx , idy + 1 ] + u[ idx , idy - 1 ] -2*u[idx,idy] )/(hy**2)
def func(xx,yy):
	return xx**3 + yy**3
def exact(xx,yy):
	return 6*xx+6*yy


#nx = int(sys.argv[1])
#ny = int(sys.argv[2])

x0 = -1 ; xn = 3
y0 = -1 ; yn = 1
'''
hy = (yn - y0)/( ny + 2 )
hx = (xn - x0)/( nx + 2 )
u = np.zeros((ny+2,nx+2))
xi = np.linspace(x0,xn,nx+2)
yi = np.linspace(y0,yn,ny+2)
xx,yy = np.meshgrid(xi,yi)

# constants
# interior indices
idx = np.arange(1,ny+1).reshape(ny,1)
idy = np.arange(1,nx+1).reshape(1,nx)
'''


nxs = [2*(2**i) for i in range(3,8)]
nys = [2**i for i in range(3,8)]
errs = []; hxs =[];hys=[]
for nx,ny in zip(nxs,nys):
	hy = (yn - y0)/( ny + 1 )
	hx = (xn - x0)/( nx + 1 )
	u = np.zeros((ny+2,nx+2))
	xi = np.linspace(x0,xn,nx+2)
	yi = np.linspace(y0,yn,ny+2)
	xx,yy = np.meshgrid(xi,yi)
	hxs.append(hx);hys.append(hy)
	idx = np.arange(1,ny+1).reshape(ny,1)
	idy = np.arange(1,nx+1).reshape(1,nx)
	lapu = lap(u,idx,idy,hx,hy)
	ue = exact(xx,yy)[idx,idy]
	errs.append(relerr(lapu,ue))
	
for i in range(1,len(errs)):
	p = np.log(errs[i-1]/errs[i])/np.log((hxs[i-1] + hys[i-1])/(hxs[i] + hys[i]))
	print(p)

#visualization
xx = xx[idx,idy]; yy = yy[idx,idy]; 
fig = plt.figure(figsize=(14,8))
ax = fig.add_subplot(1,1,1)
plt1 = ax.contourf(xx,yy,lapu,100,cmap='jet')
fig.colorbar(plt1)
plt.savefig('approx.png')

fig = plt.figure(figsize=(14,8))
ax = fig.add_subplot(1,1,1)
plt2 = ax.contourf(xx,yy,ue,100,cmap='jet')
fig.colorbar(plt2)
plt.savefig("exact.png")
