#Salman Bin Kashif
#Sarupria Group
#Clemson University
#March 13, 2018

import numpy as np
import LJ
import Init

def Calc_KE(epsilon,sigma,mass,v):
	global E_kin
	E_kin=0.0
	
	for i in range(0,Init.ip):
		v_i_sq=v[i][0]*v[i][0]+v[i][1]*v[i][1]+v[i][2]*v[i][2] #Calculating the velocity
		E_kin+=0.5*mass*v_i_sq
	return E_kin
	


