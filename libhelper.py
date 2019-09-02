import numpy as np
def relerr(uapprox,uexact):
    return np.linalg.norm((uapprox-uexact))/np.linalg.norm((uexact))

def abserr(uapprox,uexact):
    return np.linalg.norm((uapprox-uexact))

def cond(A):
    return np.linalg.cond(A)

def norm(M):
    return np.linalg.norm(M,2)

def inv(A):
    return np.linalg.inv(A)

def DM(pts,ctrs):
    return np.abs(np.subtract.outer(pts,ctrs))

def DM2D(x,y):
    return np.sqrt(np.abs(np.subtract.outer(x,x))**2 +\
                   np.abs(np.subtract.outer(y,y))**2)

def DM3D(x,y,z):
    return np.sqrt(np.abs(np.subtract.outer(x,x))**2 +\
                   np.abs(np.subtract.outer(y,y))**2 +\
                   np.abs(np.subtract.outer(z,z))**2)

def DM3D_eval(x,y,z,xe,ye,ze):
    return np.sqrt(np.abs(np.subtract.outer(xe,x))**(2)+\
                   np.abs(np.subtract.outer(ye,y))**(2)+\
                   np.abs(np.subtract.outer(ze,z))**(2))
def reshp(x):
    return x.reshape(len(x),1)

def gesol(lhs,rhs):
    return np.linalg.solve(lhs,rhs)
def reshp(x):
    return x.reshape(len(x),1)
def dm_nd(r,p):
    if len(r) == 1:
        x=r[0]
        return (np.abs(reshp(x)-reshp(x).T)**p)**(1.0/p)
    
    elif len(r) == 2:
        x=r[0].flatten();y=r[1].flatten()
        return (np.abs(reshp(x)-reshp(x).T)**p +\
        np.abs(reshp(y)-reshp(y).T)**p)**(1.0/p)
    
    elif len(r) == 3:
        x=r[0].flatten();y=r[1].flatten();z=r[2].flatten()
        return (np.abs(reshp(x)-reshp(x).T)**p +\
         np.abs(reshp(y)-reshp(y).T)**p +\
         np.abs(reshp(z)-reshp(z).T)**p)**(1.0/p)
    
    else:
        print("No Higher Dimensions Implemented")
