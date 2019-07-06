import math
i=0 
f= open("file11.txt", "w")

for i in range(2,1001):
	x = ((math.cos(i)*math.sin(i))/i)
	f.write("\n"+str(x))

#write "pebble"+str(i)
	#if i >= 15:

	#f.write("\npebble"+str(i) )
	#f.write("\n{\ntype searchableSphere;")
	
	#for i in range(0,31,30):
	#	n= str(i)
	#	f.write("\n{\ntype searchableSphere;")
	#	f.write("\ncentre("+ n +" 0 0);")
	#	f.write("\nradius 0.15\n}")


#{
 #       type searchableSphere;
  #      centre (-.15 0 0);
   #     radius 0.15;
    #}
