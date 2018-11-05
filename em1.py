import numpy as np
import math as m
import matplotlib.pyplot as plt
f= open("em1d.txt", "w")

#Forward Euler Method for ODE
#dy/dx = −y + sin(x)
#domain [0,10]

#IC
y0=2

#domain [0,10]
x0=0
xf=100

#steps
n=1000
dx=(xf-x0)/n

#y(x) & x domain
x=np.linspace(x0,xf,n)
y=np.zeros([n])
y[0]=1
#approx. solution ode
for i in range(0,n):
	y[i] = dx* (-y[i-1] + np.sin(x[i-1])) + y[i-1]

#data from ode solution
#for i in range(n):
	#print((x[i],y[i]))

	#writing data to file
 	#f.write("(" + str(x[i]) + "," + str(y[i]) + ")" + "\n" ) 
	f.write("{" + str(x[i]) + "," + str(y[i]) + "}" + "," +"\n" ) 
#x,y plot
plt.plot(x,y,'o')
plt.xlabel("X")
plt.ylabel("Y")
plt.title("Approx. ODE EulerMethod")
#plt.axis([0,10,-1,1])
#plt.show()
plt.savefig('em1.png')





