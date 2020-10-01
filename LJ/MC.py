#Salman Bin Kashif
#Sarupria Group
#Clemson University
#March 21, 2018

import numpy as np
import matplotlib.pyplot as plt
import Init
import LJ
import LJ_MC
import math
import random

def Run_MC_Metropolis(eps,sigma,mass,time,dt):
	global r #Variable for storing the position coordinates
	global r_tmp #Variable for storing the stopgap position coordinates
	global at #Variable for storing the acceleration
	global old_ener #Variable for storing inital energy
	global new_ener	#Variable for storing updated energy
	global accept # Varaible for checking whether the move was an accept or a reject
	global final_sys_ener #Variable for storing the selected energy, previus step or current step

	
	#Initializing
	Init.initialize()
	Init.print_init()

	r=Init.r

	n_atoms=Init.ip
	
	#Initializing the number of steps
	n_moves=10000

	#Providing the allowable extremes for a move
	max_translate = 0.2	

	#Providing the temperature obtained from the NVE ensemble
	T= float(input('Enter the temperature:'))

	sigma = 1.0
	eps = 1.0
	
	#Counter for number of accepts and rejects
	naccept=0
	nreject=0

	old_ener=0.0
	new_ener=0.0

	r_tmp =[0.0]*3

	#Initlializing the energy for the whole system
	old_ener=LJ.Calc_LJ(eps,sigma,r)[3]
	
	final_sys_ener=[[0.0] for i in range(0,n_moves)]
	
	Ener_sum=0.0
	Ener_avg=0.0
	
	print ('Move|Potential Energy|Average Potential Energy| Number of Accepts| Number of Rejects')
	with open ('MC_Ener.out','w') as f:
		f.write ('Move|Potential Energy|Average Potential Energy| Number of Accepts| Number of Rejects')
		for i in range(0,n_moves):
			
			#Select the random atom for the move		
			rand_atom=random.randint(0,n_atoms-1)

			#Calculate the energy only with respect to the random particle, not for whole system
			ener_rand_old=LJ_MC.Calc_LJ(eps,sigma,r,rand_atom)[3]
			
			#Storing this energy temporarily
			r_tmp[0]=r[rand_atom][0]
			r_tmp[1]=r[rand_atom][1]
			r_tmp[2]=r[rand_atom][2]

			#Determining the size of move
			del_x=random.uniform(-max_translate,max_translate)
			del_y=random.uniform(-max_translate,max_translate)
			del_z=random.uniform(-max_translate,max_translate)

			#Making the move
			r[rand_atom][0] += del_x
			r[rand_atom][1] += del_x
			r[rand_atom][2] += del_x
			
			#Calculating the updated energy only with respect to the selected particle
			ener_rand_new=LJ_MC.Calc_LJ(eps,sigma,r,rand_atom)[3]
			
			#Updating the system energy
			new_ener=old_ener-ener_rand_old+ener_rand_new

			#Initally setting the move for reject
			accept=False

			#Now checking the differnce between old and new energies and applying Metropolis criteria
			if (new_ener <= old_ener):
				accept=True 	#Making the move if crieteria satisfied
			else:
				x=math.exp(-(1/T)*(new_ener-old_ener)) 	#If energy criteria not satisfied, still accepting the move with a probability

				if(x >= random.uniform(0.0,1.0)):
					accept=True
				else:
					accept=False
			if accept==True:
				naccept+=1
				final_sys_ener[i]=new_ener #Updating the energy if move accepted
			else:
				nreject+=1

			#Setting the coordinates back to initial if move rejected
				r[rand_atom][0]=r_tmp[0]
				r[rand_atom][1]=r_tmp[1]
				r[rand_atom][2]=r_tmp[2]

				final_sys_ener[i]=old_ener
			Ener_sum+=final_sys_ener[i]
			Ener_avg=Ener_sum/(i+1) #Running average of energy
			print ('%d %20.8f %20.8f %5d %5d' % (i, final_sys_ener[i], Ener_avg ,naccept, nreject))
			f.write(('%d %20.8f %20.8f %5d %5d\n' % (i, final_sys_ener[i],Ener_avg, naccept, nreject)))
			old_ener=final_sys_ener[i]

	f.close()
#Uncomment the next line if you want to run only this file
#Run_MC_Metropolis(1.0,1.0,1.0,1,0.001)
