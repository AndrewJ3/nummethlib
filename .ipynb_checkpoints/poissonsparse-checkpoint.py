def fd2poissonsp(f,g,a,b,m):
    
    h=(b-a)/(m+1)
    
    # Sparse Differentiation Matrices (2D Laplacian)

    l1=[1]*(m-1)
    l0=[-2]*(m)

    D2=sp.diags( [l1,l0,l1], [-1,0,1], format='csr')
    
    I=sp.eye(m)
    D2y=sp.kron(I,D2)
    D2x=sp.kron(D2,I)
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
    u=spsolve(L,f)
    u=np.reshape(u,(m,m))

    return u
