def rqi(A,x0,ep):
    I=np.eye(len(A))
    lam=(x0.T @A@x0)/(x0.T @x0)
    err=norm2(A@x0-lam*x0)
    i=0
    while(err>ep):
        y=solve((A-lam*I),x0)
        x0=y/norm2(y)
        lam=(x0.T @A@x0)/(x0.T @x0)
        err=norm2(A@x0-lam*x0)
        i+=1
    print('the error converged to',err, 'in',i,'iterations')
    return x0,lam
