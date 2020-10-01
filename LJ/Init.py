#Salman Bin Kashif
#Sarupria Group
#Clemson University
#March 10, 2018

import numpy as np
import random
import math

#Generating the initial configuration

def initialize():

	global L #Box Length
	global N #Number of atoms
	global Atom_ID #Assigning an ID for a particular atom type, useful for VMD visualization
	global Atom_Type #Assigning the atom name
	global r #Variable for position coordinates
	global v #Variable for velocity coordinates 
	global V #Variable for volume
	global ip #Variable for total number of atoms inserted
	
	print ('Please specify the system. The test system was specified at 216 particles and volume as 287.496 units')
	
	N=(float)(input('Enter the number of atoms:'))
	V=(float)(input('Enter the volume:'))
	L=V**(1.0/3.0) # Calculating Box length
	print ('Box Length:%f'%(L))
	ngrids=N**(1.0/3.0) #Calculating number of grids
	print ('Number of grids:%f'%(ngrids))
	ngrids=(int)(math.floor(ngrids)+1) #Rounding number of grids to nearest integer
	print ('Integer number of grids:%f'%(ngrids))
	grid_size= L/ngrids #Calculating grid size
	print ('Grid Size:%f'%(grid_size))
	
	if grid_size<=0.5:
		print('The grid size is very small. It may lead to really large potential values and a possible atom overlap. Do you wish to continue? Enter 1 for Yes and 0 for No')
		check=input()
		if check==0:
			exit(1)
	Atom_ID=1;  #Assigning an identification to the type of item
	Atom_Type="C" # Name of the atom
	
	r=np.zeros(((int)(N),3))
	v=np.zeros(((int)(N),3))
	ip=0

	#Generation of initial positions
	print('Generating Initial coordinate')
	for i in range(0,ngrids):
		for j in range(0,ngrids):
			for k in range(0,ngrids):
				if ip<N:
					x=i*grid_size+0.5 * grid_size
					y=j*grid_size+0.5 * grid_size
					z=k*grid_size+0.5 * grid_size
					r[ip]=[x,y,z]
					print (r[ip])
					ip=ip+1
	#Seeding the initial velocity for consistent initial conditions for comparison
	random.seed(50)
	
	#Generating initial positions
	for i in range (0,ip):
		v[i]=[random.gauss(0,1) for j in range(0,3)]
		vxm = np.mean([v[i][0] for i in range(0,ip)])
		vym = np.mean([v[i][1] for i in range(0,ip)])
		vzm = np.mean([v[i][2] for i in range(0,ip)])
		for i in range(0,ip):
			v[i][0] = v[i][0] - vxm
			v[i][1] = v[i][1] - vym
			v[i][2] = v[i][2] - vzm
	print ('Number of atoms finally in grid space:%d'%(ip))
	return r,v

#writing the positions and velocities to file
def print_init():
	with open('Init_Cond.xyz','w') as f:
		f.write("%d\nInitial Conditions\n" % (N))
		for i in range(0,ip):
			f.write('%s%8.3f%8.3f%8.3f\n' % (Atom_Type,r[i][0], r[i][1], r[i][2]))
			
	with open('Init_Vel.out','w') as e:
		for i in range(0,ip):
			e.write('%10.6f%10.6f%10.6f\n' % (v[i][0], v[i][1], v[i][2]))
			print('%10.6f%10.6f%10.6f\n' % (v[i][0], v[i][1], v[i][2]))
	e.close()
	f.close()

