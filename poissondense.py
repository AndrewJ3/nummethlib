def fd2poisson(f,g,a,b,m):
    
    h=(b-a)/(m+1)
    
    # Dense Differentiation Matrices (2D Laplacian)

    z1=np.array([-2,1])
    z2=np.zeros(m-2)
    z=np.block([z1,z2])
    D2=la.toeplitz(z,z)
    I=np.eye(m)
    D2x=np.kron(D2,I)
    D2y=np.kron(I,D2)
    L=(h**-2)*(D2x+D2y)
    
    #Boundary Conditions
    ubs=g[0,1:m+1]
    ubn=g[m+1,1:m+1]
    ube=g[1:m+1,m+1]
    ubw=g[1:m+1,0]

    f[:,0]-=ubw*(h**-2)
    f[:,m-1]-=ube*(h**-2)
    f[0,:m]-=ubs*(h**-2)
    f[m-1,:m]-=ubn*(h**-2)

    
    f=np.reshape(f,m*m)
    
    #Solving the System
    u=np.linalg.solve(L,f)
    u=np.reshape(u,(m,m))

    return u
