i=0
n=0 
f= open("sphere.txt", "w")
while 1==1:
	i=i+15

	#write "pebble"+str(i)
	n=n+1

	if n >=2:
		if i >= 3:
			break
	x1=str(0)
	x2=str(3.6)
	x3=str(0)
	x4=str(6)
	x5=str(1.2)
	x6=str(4.2)
	x7=str(4.8)
	f.write("\npebble"+str(n))
	f.write("\n{\n	type searchableSphere;")
	f.write("\n	centre ("+x1 +" "+ x2 +" "+ x3 +");")
	f.write("\n	radius 0.3	\n}")

	f.write("\npebble"+str(n))
	f.write("\n{\n	type searchableSphere;")
	f.write("\n	centre ("+x1 +" "+ x6 +" "+ x3 +");")
	f.write("\n	radius 0.3	\n}")

	f.write("\npebble"+str(n))
	f.write("\n{\n	type searchableSphere;")
	f.write("\n	centre ("+x1 +" "+ x7 +" "+ x3 +");")
	f.write("\n	radius 0.3	\n}")

	f.write("\npebble"+str(n))
	f.write("\n{\n	type searchableSphere;")
	f.write("\n	centre ("+x1 +" "+ x2 +" -."+ x4 +");")
	f.write("\n	radius 0.3	\n}")

	f.write("\npebble"+str(n))
	f.write("\n{\n	type searchableSphere;")
	f.write("\n	centre ("+x1 +" "+ x6 +" -."+ x4 +");")
	f.write("\n	radius 0.3	\n}")

	f.write("\npebble"+str(n))
	f.write("\n{\n	type searchableSphere;")
	f.write("\n	centre ("+x1 +" "+ x7 +" -."+ x4 +");")
	f.write("\n	radius 0.3	\n}")

	f.write("\npebble"+str(n))
	f.write("\n{\n	type searchableSphere;")
	f.write("\n	centre ("+x1 +" "+ x2 +" -"+ x5 +");")
	f.write("\n	radius 0.3	\n}")

	f.write("\npebble"+str(n))
	f.write("\n{\n	type searchableSphere;")
	f.write("\n	centre ("+x1 +" "+ x6 +" -"+ x5 +");")
	f.write("\n	radius 0.3	\n}")

	f.write("\npebble"+str(n))
	f.write("\n{\n	type searchableSphere;")
	f.write("\n	centre ("+x1 +" "+ x7 +" -"+ x5 +");")
	f.write("\n	radius 0.3	\n}")



	if i >= 3:
		f.close
	
