def house(A):
    
    #copy A into R to avoid overwriting A
    m,n=A.shape
    R=np.zeros((m,n))
    W=np.zeros((m,n))
    for i in range(0,m):
        R[i]=A[i]

    #householder transformation
    for i in range(0,n):
        x=R[i:m+1,i]
        e=np.zeros(len(x))
        e[0]=1.0        
        u=np.sign(x[0])*norm2(x)*e+x
        u/=norm2(u)

        # finding R and W
        u=np.array([u]).T
        uT=u.T
        R[i:m+1,i:n+1] -= 2.0*u@(uT@R[i:m+1,i:n+1])
        W[i:m,i]=u[:,0]
        
    return R[:n,:n],W
