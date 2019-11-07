def input_args():	
	file = open("input.in","r")
	print(file.readline())
	print(file.readline())
	print(file.readline())
	global nx,ny,finaltime, dt,tol
	for line in file:
		word,num=line.split(" ")
		if word=="nx":
			nx = int(num)
		elif word=="ny":
			ny = int(num)
		elif word=="finaltime":
			finaltime = float(num)
		elif word=="dt":
			dt = float(num)
		else:
			tol= float(num)
		print(word+" "+num)

input_args()
print(nx)
