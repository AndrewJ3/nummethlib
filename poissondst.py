def fd2poissondst(f,g,a,b,m):
    
    h=(b-a)/(m+1)
    
    #Boundary Conditions
    ubs=g[0,1:m+1]
    ubn=g[m+1,1:m+1]
    ube=g[1:m+1,m+1]
    ubw=g[1:m+1,0]

    f[:m+1,0]-=ubw*(h**-2)
    f[:m+1,m-1]-=ube*(h**-2)
    f[0,:m+1]-=ubs*(h**-2)
    f[m-1,:m+1]-=ubn*(h**-2)

    fhat=dst(dst(f,1).T,1)/(2*m+2)
    
    j=np.arange(1,m+1)
    
    denom=np.add.outer(np.cos(j*np.pi/(m+1)),np.cos(j*np.pi/(m+1))) -2
    uhat=((h**2)*fhat)/(2*denom)
    
    u=idst(idst(uhat,1).T,1)/(2*m+2)

    return u
