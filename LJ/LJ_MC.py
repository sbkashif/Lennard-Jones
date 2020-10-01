#Salman Bin Kashif
#Sarupria Group
#Clemson University
#March 13, 2018

import numpy as np
import random
import Init
import math

#For Monte Carlo simulation
#The code is similar to LJ file with the only difference being calculation of energy only with respect to a particle
#Follow the comments in LJ.py file for details

def Calc_LJ(eps,sigma,r,j):	
	global sigma2
	global fx
	global fy
	global fz
	global E_pot
	global rc
	
	sigma2=sigma**2

	fx=[0.0]*len(r)
	fy=[0.0]*len(r)
	fz=[0.0]*len(r)
	E_pot=0.0

	for i in range(0,Init.ip):
		if i!=j: #Just avoid the same particle as we do not know the index beforehand
			dx=r[i][0]-r[j][0]
			dy=r[i][1]-r[j][1]
			dz=r[i][2]-r[j][2]
			
			dx=dx - round(dx/Init.L)*Init.L
			dy= dy -round(dy/Init.L)*Init.L
			dz= dz -round(dz/Init.L)*Init.L

									
			r2=dx**2+dy**2+dz**2

			rc= 5.0 * sigma

			if (r2 == 0):
				print ("Error dist=0")
				print ('i=%d j=%d,xi=%f yi=%f zi=%f xj=%f yj= %f zj=%f' % (i,j,r[i][0],r[i][1],r[i][2],r[j][0],r[j][1],r[j][2]))
				print ('dx=%f, dy=%f, dz=%f, r2=%f' % (dx,dy,dz,r2))
				print ('%f' % (Init.L))
				print('Particle overlap reported')
				exit(1)
			elif r2 < (rc):
				dist=math.sqrt(r2)
				fpr=4.0 * eps * (12.0*(sigma**12)*(dist**(-13))-6.0*(sigma**6)*(dist**(-7)))/dist
				fxi=fpr * dx
				fyi=fpr * dy
				fzi=fpr * dz
				E_pot_i=4.0*eps*((sigma **12)/(dist**12)-(sigma**6)/(dist **6))
				E_pot_rc=Pot_Rc(rc, eps, sigma)
				E_pot_i=E_pot_i-E_pot_rc
			else:
				fxi=0
				fyi=0
				fzi=0
				E_pot_i=0
				
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

def Pot_Rc(rc, eps, sigma):
	dis=math.sqrt(rc)
	E_pot_rc=4.0 *eps*((sigma**12)/(dis**12)-(sigma**6)/(dis**6))
	return E_pot_rc
