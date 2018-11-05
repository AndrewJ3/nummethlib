#Forward Euler Method for ODE
#dy/dx = −y + cos(x)
#domain [0,10]

import numpy as np
import math as m
import matplotlib.pyplot as plt

#f= open("em1de.txt", "w")
f1 = open("emerror.txt", "w")

#IC
y0=1

#domain [0,10]
x0=0
xf=100

#steps
n=1000
dx=(xf-x0)/n

#y(x) & x domain
x=np.linspace(x0,xf,n)
y=np.zeros([n])
yte=np.zeros([n])
ye=np.zeros([n])
#y[0]=y0
#approx. solution ode
for i in range(0,n):
	y[0]=y0
	y[i] = dx* (-y[i-1] + np.cos(x[i-1])) + y[i-1]
	
#truncation error 2nd order
	yte[i]=(1/(i+1))*(dx**(i))*(-y[i-1]+np.cos(x[i-1]))
#error  em1we2.py
	#ye[i]=y[i]+yte[i]

#analytic solution
#for i in range(0,n):
	ya=0.5*np.exp(-x)*(1+np.exp(x)*np.cos(x)+np.exp(x)*np.sin(x))
#1/2 E^-x (1 + E^x Cos[x] + E^x Sin[x])
	#Error
	#x=0.3
	#ye = np.abs(ya-y[i])
	#print(ye)
	#print(i)
  #data from ode solution
  #for i in range(n):
	#print((x[i],y[i]))
	#writing data to file
 	#f.write("(" + str(x[i]) + "," + str(y[i]) + ")" + "\n" ) 

	#f1.write("(" + str(dx) + "," + str(ye) + ")" + "\n" ) 

#x,y plot
plt.plot(x,ye,'o-')
plt.plot(x,ya,'o-')
plt.xlabel("X")
plt.ylabel("Y")
	#plt.title("Approx. ODE EulerMethod")
	#plt.axis([0,10,-1,1])
	#plt.show()
plt.savefig('emnvawtr2.png')
 
#error plot
#plt.plot(i,ye,'o')
	#plt.axis([0,100,0,1])
#plt.savefig('em1error.png')


