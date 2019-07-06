import math
x=3
i=1
f= open("Geo.txt", "w")
while 1==1:
	i=i+1
	if i> 50:
		#z = (math.pow(x,i))
		break
	z= (math.pow(x,i))
	f.write("\n"+ str(z))

#f= open("file11.txt", "w")

#for i in range(2,1001):
#	x = ((math.cos(i)*math.sin(i))/i)
#	f.write("\n"+str(x))

