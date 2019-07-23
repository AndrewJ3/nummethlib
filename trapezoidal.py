def func(x):
    return np.cos(x)

def gfunc(x):
    return np.sin(x)

def integral_func(g,b,a):
    return g(b) - g(a)

a=0;b=1
errs=[]
ns=[100*i for i in range(1,10)]
for n in ns:
    x=np.linspace(a,b,n)
    u=np.zeros(n)
    dx=(b-a)/n
    u=dx*((func(a)+func(b))/2 + np.sum(func(x[1:-1])))
    errs.append(np.abs(u-integral_func(gfunc,b,a)))
    
# plt.loglog(ns,errs)
# plt.plot(x,np.sin(x))
area_under_curve = np.zeros((n**2,2))
for i in range(n**2):
    ptx0 = (b)*np.random.rand(1) + a
    pty0 = (gfunc(x[-1]) - gfunc(x[0]))*np.random.rand(1) + gfunc(x[0])
    ii = np.random.randint(0,n)
    if ptx0 >= x[ii] and pty0 <= gfunc(x[ii]):
        area_under_curve[ii,:] = np.array([ptx0,pty0]).flatten()
