import math
#i=0 
#def sum (arg1, arg2):
	#f= open("file12.txt", "w")
	#total = arg1 + arg2; # Here total is local variable.
	#x=sum(1,2)
	#print(x)
	#f.write("\n"+str(x))
f= open("file12.txt", "w")
sum = 0
for n in range(1,100):
   sum = (math.cos(n)*math.sin(n))/(n)
f.write("\n"+str(sum))

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
#math.cos(i)*math.sin(i))/i)
