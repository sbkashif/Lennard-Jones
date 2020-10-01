#Salman Bin Kashif
#Sarupria Group
#Clemson University
#March 13, 2018

import numpy as np
import random
import Init
import math

def dW(delta_t):
	return np.random.normal(loc=0.0,scale=np.sqrt(delta_t))


#Calculating L-J potential
def Calc_LJ(eps,sigma,r,v,sig,gamma,del_t):	
	global sigma2
	global fx #x-coordinate for force
	global fy #y-coordinate for force
	global fz #z-coordinate for force
	global E_pot #Variable to store potential energy
	global rc #Variable for cut-off radius
	
	sigma2=sigma*sigma

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
			dy=dy -round(dy/Init.L)*Init.L
			dz=dz -round(dz/Init.L)*Init.L

			vijx=v[i][0]-v[j][0]
			vijy=v[i][1]-v[j][1]
			vijz=v[i][2]-v[j][2]
			dot_prod=vijx*dx+vijy*dy+vijz*dz			

						
			r2=dx*dx+dy*dy+dz*dz

			rc= 2.5 * sigma
			sqrt_dt=math.sqrt(del_t)
			
			#Termination step if particle overlap is encountered as it is not possible in real system

			if (r2 == 0):
				print ("Error dist=0")
				print ('i=%d j=%d,xi=%f yi=%f zi=%f xj=%f yj= %f zj=%f' % (i,j,r[i][0],r[i][1],r[i][2],r[j][0],r[j][1],r[j][2]))
				print ('dx=%f, dy=%f, dz=%f, r2=%f' % (dx,dy,dz,r2))
				print('Particle overlap reported')
				exit(1)
			elif r2 < (rc):
				#print ('Inside cut-off')
				dist=math.sqrt(r2) #Inter-particle distance
				fpr=4 * eps * (12.0*(sigma**12)*(dist**(-13))-6.0*(sigma**6)*(dist**(-7)))/dist #Calculating forces
				fxi=(fpr * dx)
				#fx[j]=(fpr * dx)/r2
				fyi=(fpr * dy)
				#fy[j]-=(fpr * dy)/r2
				fzi=(fpr * dz)
				#fz[j]-=(fpr * dz)/r2
				E_pot_i=4.0*eps*((sigma **12)/(dist**12)-(sigma**6)/(dist **6))
				E_pot_rc=Pot_Rc(rc, eps, sigma) #Calculating potential shift distance
				E_pot_i=E_pot_i-E_pot_rc #Applying potential shift truncationi

#Calculating the random forces				

				fd=-gamma*dot_prod/r2
				fxi1=(dx*fd)
				#fx[j]-=(dx*fd)
				fyi1=(dy*fd)
				#fy[j]-=(dy*fd)
				fzi1=(dz*fd)
				#fz[j]-=(dz*fd)
				
				fr=sig/dist/sqrt_dt
				z=random.gauss(0.0,1.0)*fr
				fxi2=(dx*z)
				#fx[j]-=(dx*z)
				#fr=sig/dist/sqrt_dt
				z=random.gauss(0.0,1.0)*fr
				fyi2=(dy*z)
				#fy[j]-=(dy*z)
				#fr=sig/dist/sqrt_dt
				z=random.gauss(0.0,1.0)*fr
				fzi2=(dz*z)
				#fz[j]-=(dz*z)
				
			else:
				#Condition for particle outside cut-off
				#print('Outside cutoff')
				fxi=0
				fxi1=0
				fxi2=0
				fyi=0
				fyi1=0
				fyi2=0
				fzi=0
				fzi1=0
				fzi2=0
				E_pot_i=0

			#Updating the forces and the potential energy
			fx[i]=fx[i]+fxi+fxi1+fxi2
			fx[j]=fx[j]-fxi-fxi1-fxi2
			fy[i]=fy[i]+fyi+fyi1+fyi2
			fy[j]=fy[j]-fyi-fyi1-fyi2
			fz[i]=fz[i]+fzi+fzi1+fzi2
			fz[j]=fz[j]-fzi-fzi1-fzi2
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
