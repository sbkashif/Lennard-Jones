#Salman Bin Kashif
#Sarupria Group
#Clemson University
#March 22, 2018

#Importing all the modules

import Init
import LJ
import Acceleration
import Vel_verlet_NVT
import Vel_verlet_NVE
import Leap_Frog_NVT
import Leap_Frog_NVE
import MC
import numpy
import random


def run():
	
	#Selecting the type of simulation
	take_0=input('Please select the type of run | 1- MD 2- MC:')
	if take_0==1:
	#Selecting the integrator	
		take=input('Please select the integrator | 1- Velocity Verlet 2- Leap Frog :')
	#Selecting the ensemble
		if take==1:
			take_2=input('Please select the ensemble | 1- NVE 2- NVT :')
			if take_2==1:
				Vel_verlet_NVE.Run_Vel_ver(1.0,1.0,1.0,1.0,0.001)
			if take_2==2:
				take_3=input('Please input the average temperature you got from NVE or put a temperature of your choice(in dimensionless units as per L-J phase diagram:')
				t_nvt_vv=(float)(take_3)
				Vel_verlet_NVT.Run_Vel_ver(1.0,1.0,1.0,1.0,0.001,0.001,t_nvt_vv)
			else:
				print('Please select a valid option. Run the program again')
				exit(1)
		if take==2:
			take_4=input('Please select the ensemble | 1- NVE 2- NVT :')
			if take_4==1:
				Leap_Frog_NVE.Run_Leap_Frog(1.0,1.0,1.0,1.0,0.001)
			elif take_4==2:
				take_5=input('Please input the average temperature you got from NVE or put a temperature of your choice(in dimensionless units as per L-J phase diagram:')
				t_nvt_lf=(float)(take_5)
				Leap_Frog_NVT.Run_Leap_Frog(1.0,1.0,1.0,1.0,0.001,t_nvt_lf,0.001)
			else:
				print('Please select a valid option. Run the program again')
				exit(1)
		else:
			print('Please select a valid option. Run the program again')
			exit(1)
	elif take_0==2:
		MC.Run_MC_Metropolis(1.0,1.0,1.0,1.0,0.001)
	else:
		print('Please select a valid option. Run the program again')
		exit(1)
run()
