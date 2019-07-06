i=0 
f= open("file9.txt", "w")

for i in range(0,33):
	#write "pebble"+str(i)
	#if i >= 15:
	i=i+1
	f.write("\npebble"+str(i) )
	f.write("\n{\ntype searchableSphere;")
	break
for s in range(0,31,15):
		n= str(i)
		f.write("\ncentre (."+ n +" ."+ n +" ."+ n +" "+");")
		f.write("\nradius 0.15\n}")


#{
 #       type searchableSphere;
  #      centre (-.15 0 0);
   #     radius 0.15;
    #}
