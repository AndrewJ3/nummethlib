i=0
n=0 
f= open("sph.txt", "w")
while 1==1:
	
	n=n+1

	if n >=82:
		
		break
	f.write("\npebble"+str(n))
	f.write("\n{\n	type  fixedValue;")
	f.write("\n	value uniform 1800;}")
	#f.write("\n	        Prt           0.85;")
	#f.write("\n	       value           uniform 1600;	\n}")		
	
