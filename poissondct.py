import numpy as np
def fd2poissondct(f,a,b,m):
    
    h=(b-a)/(m+1)
    
    fhat=dct(dct(f,1).T,1)/(2*m+2)
    j=np.arange(0,m+2)
    uhat = np.zeros((m,m))
    dnm=np.cos((j*np.pi)/(m+1))+np.cos((j*np.pi)/(m+1))-2

    dnm[0]=2
    uhat=(fhat*h**(2))/(2*dnm)
    uhat[0,0]=0

    u=idct(idct(uhat,1).T,1)/(2*m+2)
    
    return u,fhat,denom
