
import numpy as np
import matplotlib.pyplot as plt
import sys

npzfile = np.load("test.npz")
uvpxy=[]
for i in range(len(npzfile.files)):
	uvpxy.append(npzfile[npzfile.files[i]])
u,v,p,x,y=uvpxy

dp = 1
h =1
l= 1
ue = 0.5*(1/0.1)*((dp)/l) *(h*y[:,-1] - y[:,-1]**2)

plt.plot(y[:,-1],u[:,-1],'bv')
plt.plot(y[:,-1],ue)
plt.savefig("flow.png")
