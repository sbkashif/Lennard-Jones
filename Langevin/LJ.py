#Salman Bin Kashif
#Sarupria Group
#Clemson University
#March 13, 2018

import numpy as np
import random
import Init
import math
#Calculating L-J potential
def Calc_LJ(eps,sigma,r):	
	global sigma2
	global fx #x-coordinate for force
	global fy #y-coordinate for force
	global fz #z-coordinate for force
	global E_pot #Variable to store potential energy
	global rc #Variable for cut-off radius
	
	sigma2=sigma**2

	#Initializing force and potential energy variables

	fx=[0.0]*len(r) 
	fy=[0.0]*len(r)
	fz=[0.0]*len(r)
	E_pot=0.0

	for i in range(0,Init.ip):
		for j in range(i+1,Init.ip): #Condition to avoid same atom distnce calculation
			dx=r[i][0]-r[j][0]
			dy=r[i][1]-r[j][1]
			dz=r[i][2]-r[j][2]
			
			#Applying periodic boundary conditions
			dx=dx - round(dx/Init.L)*Init.L
			dy= dy -round(dy/Init.L)*Init.L
			dz= dz -round(dz/Init.L)*Init.L

									
			r2=dx**2+dy**2+dz**2

			rc= 5.0 * sigma
			
			#Termination step if particle overlap is encountered as it is not possible in real system

			if (r2 == 0):
				print ("Error dist=0")
				print ('i=%d j=%d,xi=%f yi=%f zi=%f xj=%f yj= %f zj=%f' % (i,j,r[i][0],r[i][1],r[i][2],r[j][0],r[j][1],r[j][2]))
				print ('dx=%f, dy=%f, dz=%f, r2=%f' % (dx,dy,dz,r2))
				print('Particle overlap reported')
				exit(1)
			elif r2 < (rc):
				dist=math.sqrt(r2) #Inter-particle distance
				fpr=4.0 * eps * (12.0*(sigma**12)*(dist**(-13))-6.0*(sigma**6)*(dist**(-7)))/dist #Calculating forces
				fxi=fpr * dx
				fyi=fpr * dy
				fzi=fpr * dz
				E_pot_i=4.0*eps*((sigma **12)/(dist**12)-(sigma**6)/(dist **6)) #Calculating potential energy
				E_pot_rc=Pot_Rc(rc, eps, sigma) #Calculating potential shift distance
				E_pot_i=E_pot_i-E_pot_rc #Applying potential shift truncation
			else:
				#Condition for particle outside cut-off
				fxi=0
				fyi=0
				fzi=0
				E_pot_i=0

			#Updating the forces and the potential energy
			fx[i]=fx[i]+fxi
			fx[j]=fx[j]-fxi
			fy[i]=fy[i]+fyi
			fy[j]=fy[j]-fyi
			fz[i]=fz[i]+fzi
			fz[j]=fz[j]-fzi
			E_pot=E_pot+E_pot_i

		with open('Interact_Force.out','w') as f:
			f.write('Lennard Jones potential calculation\n')	
			for i in range (0,len(r)):
				f.write('{:<20.6f},{:<20.6f},{:<20.6f}\n'.format(fx[i],fy[i],fz[i]))
		f.close()

	return fx,fy,fz,E_pot

#Function for calculating potential shift correction
def Pot_Rc(rc, eps, sigma):
	dis=math.sqrt(rc)
	E_pot_rc=4.0 *eps*((sigma**12)/(dis**12)-(sigma**6)/(dis**6))
	return E_pot_rc
