def vandmn(x,n):
    A=np.ones((len(x),n))
    for i in range(1,n):
        A[:,i]=x**(i)
    return A
