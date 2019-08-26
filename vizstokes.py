import numpy as np
import matplotlib.pyplot as plt
import sys

npzfile = np.load("test.npz")
uvpxy=[]
for i in range(len(npzfile.files)):
	uvpxy.append(npzfile[npzfile.files[i]])
u,v,p,x,y=uvpxy
plt.plot(x[0,:],u[0,:])
plt.savefig("flow.png")
