def house2q(W):
    m,n=W.shape
    Q=np.zeros((m,m))
    I=np.eye(m)
    for i in range(m):
        b=I[:,i]
        for j in range(n):
            w=W[:,j].reshape(m,1)
            b=b.reshape(m,1)
            gamma= -2.0*(w.T@b)
            b=b+gamma*w
        Q[i,:]=b.T
    return Q[:m,:n]
