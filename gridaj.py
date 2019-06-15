def gridaj(a,b,nx,ny):
    hx=(b-a)/(nx-1);hy=(b-a)/(ny-1)
    x=np.zeros(nx)
    y=np.zeros(ny)
    gridx=np.zeros((nx)**2)
    gridy=np.zeros((ny)**2)
    for i in range(0,nx):
        x[i] = a + i*hx
        y[i] = a + i*hy
    for i in range(0,nx):  
        for j in range(0,nx):
            gridx[i*(nx)+j] = x[j]
    
    for i in range(0,ny):  
        for j in range(0,ny):
            gridy[j*(ny)+i] = y[j]

    return gridx,gridy
