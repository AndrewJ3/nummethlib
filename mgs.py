def mgs(A):
    m,n=A.shape
    Q=np.zeros((m,n))
    R=np.zeros((n,n))
    for k in range(n):
        Q[:,k]=A[:,k]
        for i in range(k):
            R[i,k] = Q[:,i].reshape(m,1).T @Q[:,k]
            Q[:,k] -= R[i,k]*Q[:,i]
        R[k,k]=norm2(Q[:,k])
        Q[:,k]/=R[k,k]
    return Q,R
