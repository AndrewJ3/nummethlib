i=0 
f= open("schd3.txt", "w")

for i in range(0,13,1):
	#write "pebble"+str(i)
	#if i >= 15:

	#f.write("\npebble"+str(i) )
	#f.write("\n{\ntype searchableSphere;")
	
	for s in range(0,60,15):
		h= str(i)
		m= str(s)
		f.write("\n"+h+ ":"+m+"pm")

#f.write("\ncentre("+ n +" 0 0);")
		#f.write("\nradius 0.15\n}")


#{
 #       type searchableSphere;
  #      centre (-.15 0 0);
   #     radius 0.15;
    #}
