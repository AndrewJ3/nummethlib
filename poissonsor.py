#Poisson Solver using Succesive Over Relaxation
#omega=w
# file=open("prob1.txt","w")
def fd2poissonsor(f,g,a,b,m):
    tol=1e-8
    h=(b-a)/(m+1)
    wopt=2/(1+np.sin(np.pi*h))
#     print('Optimal Omega',wopt)
    w=wopt
    u=np.zeros((m+2,m+2))
    res2=np.zeros((m+2,m+2))
    ubs=g[0,1:m+2]
    ubn=g[m,1:m+2]
    ube=g[1:m+2,m]
    ubw=g[1:m+2,0]
    u[0,1:m+2]=ubs
    u[m,1:m+2]=ubn
    u[1:m+2,m]=ube
    u[1:m+2,0]=ubw
    itermax=m**2  
    for i in range(itermax):
        res=0
        for j in range(1,m):
            for k in range(1,m):
                un=(1-w)*u[j,k]+w*.25*(u[j+1,k]+u[j-1,k]+u[j,k+1]+u[j,k-1]-(h**(2))*f[j,k])
                u[j,k]=un
        for j in range(1,m):
            for k in range(1,m):
                res=res+(-4*u[j,k]+(u[j+1,k]+u[j-1,k]+u[j,k+1]+u[j,k-1]-h**(2)*f[j,k]))**2
                res2[j,k] = res
        res=np.sqrt(res2[j,k])
        if res < tol:
#             print(res)
#             print('Finished after %g Iterations'%i)
#             print('Final Residual= %g'%res)
            break;   

    return u
