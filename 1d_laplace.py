n=100;a=-1;b=1
alpha=3;beta=2
x=np.linspace(a,b,n)

# analytical solution u=(beta-alpha)x-alpha
bcrhs=np.array([[a,1],
          [b,1]])
bclhs=np.array([[alpha],[beta]])
cst=gesol(bcrhs,bclhs)

u_exact=cst[0]*x+cst[1]

# second order finite difference method iterative method SOR and Jacobi
un=np.zeros(n)
res=np.zeros(n)
h=(b-a)/(n-1)
f=np.zeros(n)
un[0]=alpha;un[n-1]=beta;
tol=1e-4
omega=1.961
w1=1.98;w0=1.93
w=(w1-w0)*np.random.random(100)+w0
iters=0
for ww in w:
    for j in range(n**3):
        for i in range(1,n-1):
            un[i] = (1-ww)*un[i] + ww*((un[i-1]+un[i+1])/2.0 - ((h**2)*f[i])/2.0)

        for i in range(1,n-1):
            res[i] = un[i] - ((1-ww)*un[i] + ww*((un[i-1]+un[i+1])/2.0 - ((h**2)*f[i])/2.0))

        if norm(res)<tol or j>1000:
#             print("iterations = ",j)
#             print("||r||_2 = ",norm(res))
            break;
    iters+=1
    if norm(res)<1e-5:
        print("found omega",ww,"in %d iterations"%iters)
        break;
        
for j in range(n**3):
        for i in range(1,n-1):
            un[i] = (1-ww)*un[i] + ww*((un[i-1]+un[i+1])/2.0 - ((h**2)*f[i])/2.0)

        for i in range(1,n-1):
            res[i] = un[i] - ((1-ww)*un[i] + ww*((un[i-1]+un[i+1])/2.0 - ((h**2)*f[i])/2.0))

        if norm(res)<tol or j>1000:
            print("iterations = ",j)
            print("||r||_2 = ",norm(res))
            break;
# for j in range(n**3):
#     for i in range(1,n-1):
#         un[i]=(un[i-1]+un[i+1])/2.0 - ((h**2)*f[i])/2.0
        
#     for i in range(2,n-1):
#         res[i]=un[i]-((un[i-1]+un[i+1])/2.0 - (((h**2)*f[i])/2.0))
    
#     if norm(res)<tol or j>10000:
#         print(j)
#         break;
        

# Absolute and Relative L2 Error
# print(abs(un-u_exact),"\n")

print("L2 Error",relerr(u_exact,un))
