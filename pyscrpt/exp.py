import numpy as np
x=1
i=0
#f= open("exp.txt", "w")
while 1==1:
	i=i+1
	if i> 100:
		#z = (math.pow(x,i))
		break
	z=(x**i) / (i)
	zn=np.sum(z)
print(zn)
	#f.write("\n"+ str(z))
	

