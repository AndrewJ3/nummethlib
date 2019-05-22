def clgs(A):
    m,n=A.shape
    Q=np.zeros((m,n))
    R=np.zeros((n,n))
    for k in range(n):
        Q[:,k]=A[:,k]
        if k != 0:
            J=np.arange(0,k)
            Q[:,k-1]=np.array([Q[:,k-1]])
            R[J,k] = Q[:,k-1].T @Q[:,k]
            Q[:,k] -= Q[:,J]@R[J,k]
        R[k,k] = norm2(Q[:,k])
        Q[:,k] /= R[k,k]
    return Q,R
