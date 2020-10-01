#Salman Bin Kashif
#Sarupria Group
#Clemson University
#March 15,2018

import numpy as np
import LJ, Init

def Calc_a(eps,sigma, mass,r,v,sig,gamma,dt):
	global ax
	global ay
	global az
	
	Fx,Fy,Fz,E_p=LJ.Calc_LJ(eps,sigma,r,v,sig,gamma,dt)
	a=[[0.0]*3 for i in range(0,Init.ip)]
	for i in range(0,Init.ip):
		a[i][0]=Fx[i]/mass
		a[i][1]=Fy[i]/mass
		a[i][2]=Fz[i]/mass
	with open('Acceleration.out','w') as f:
		f.write('Acclearation of all the particles at given time step\n')
		for i in range(0, Init.ip):
			f.write('{:<20.6f},{:<20.6f},{:<20.6f}\n'.format(a[i][0],a[i][1],a[i][2]))
	f.close()
	return a

	
				
