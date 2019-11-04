import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import sympy as sym
import sys

n = int(sys.argv[1])
dim = int(sys.argv[2])
tol = 1e-14

#=========================(==============================================
def ufoo(x,y,z,pt,d):
	px,py,pz = pt
	return 1/((x-px+1e-14)**2+\
			(y-py+1e-14)**2+\
			(z-pz+1e-14)**2)**((d-2)/2)

def dufoo():
	x,y,z,xj,yj,zj,d = sym.symbols('x y z xj yj zj d')
	pe = 1.0/((x-xj)**2+(y-yj)**2+(z-zj)**2)**((d-2)/2)
	df = sym.Matrix([sym.diff(pe,x),sym.diff(pe,y),sym.diff(pe,z)])
	return sym.lambdify((x,y,z,xj,yj,zj,d),sym.simplify(df[0]),'numpy'),\
		sym.lambdify((x,y,z,xj,yj,zj,d),sym.simplify(df[1]),'numpy'),\
		sym.lambdify((x,y,z,xj,yj,zj,d),sym.simplify(df[2]),'numpy')


def newton(x0,y0,z0,pt,f,df,d,maxiter,epsilon):
	xk = x0;yk = y0;zk = z0;
	for k in range(maxiter):
		fxk = f(xk,yk,zk,pt,d)
		dfxk = df[0](xk,yk,zk,pt[0],pt[1],pt[2],d)
		dfyk = df[1](xk,yk,zk,pt[0],pt[1],pt[2],d)
		dfzk = df[2](xk,yk,zk,pt[0],pt[1],pt[2],d)
		if abs(fxk) < epsilon:
#			print('Found solution after',k,'iterations.')
			nrm = np.sqrt(xk**2+yk**2+zk**2)
			return xk/nrm,yk/nrm,zk/nrm
		elif (dfxk == 1e-14):
#			print('Zero derivative. No solution found.')
			return None
		xk = xk - fxk/dfxk
		yk = yk - fxk/dfyk
		zk = zk - fxk/dfzk
#	print('Exceeded maximum iterations. No solution found.')
#=========================(==============================================

th = (2*np.pi)*np.random.random(n) 
lam = (np.pi)*np.random.random(n)
pts = np.zeros((n**2,3))

pts[:,0] = np.outer(np.cos(th),\
        np.sin(lam)).flatten()
pts[:,1] = np.outer(np.sin(th),\
             np.sin(lam)).flatten()
pts[:,2] = np.outer(np.ones(n),\
            np.cos(lam)).flatten()

minpts = np.zeros((n**2,3))
for i in range(n**2-1):
	x,y,z=newton(pts[i,0],\
				pts[i,1],\
				pts[i,2],\
				pts[i+1,:],\
                 ufoo,\
                 dufoo(),dim,500,tol)

	minpts[i,0] = x
	minpts[i,1] = y
	minpts[i,2] = z
print(minpts.shape)

fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')
ax.scatter(minpts[:,0],minpts[:,1],minpts[:,2])
plt.savefig('fig.png')
