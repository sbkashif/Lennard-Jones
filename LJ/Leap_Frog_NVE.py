#Salman Bin Kashif
#Sarupria Group
#Clemson University
#March 20, 2018

import numpy as np
import Init
import Acceleration
import LJ
import Kinetic_Energy
import pylab

#Leap Frog Algorithm for NVE Thermostat

def Run_Leap_Frog(eps,sigma,mass,time,dt):
	global rt
	global vt
	global at
	global E_potential
	global E_kinetic

	#Initialization
	Init.initialize()
	Init.print_init()

	rt=Init.r #Assigning inital values to positions
	vt=Init.v #Assigning initial values to velocity
	at=[[0.0]*3 for i in range(0,Init.ip)] #Initializing the accelaration
	sum_a=[[0.0]*3 for i in range(0,Init.ip)] #Initilzing the a(t) and a(t+dt) varialble
	dth=dt*0.5 #Half time step
	nsteps=(int)(time/dt) #Number of steps to run the simulation
	#Initialization of energy and temperature variables
	E_potential=[0.0]*nsteps
	E_kinetic=[0.0]*nsteps
	E_total=[0.0]*nsteps
	T_avg=0.0
	T_sum=0.0

	at=Acceleration.Calc_a(eps,sigma,mass,rt) #Assigning initial values to acceleration by calling the function
	
	E_kin=[0.0]*(nsteps)
	T=[0.0]*(nsteps)
	print('Time step | Potential Energy | Kinetic Energy| Total Energy|Temperature|Average Temperature\n')
	with open('Leap_Frog.out','w') as f, open('Leap_Frog_Energy','w') as e, open('VMD_pos_leapfrog.xyz','w') as d: #Initializing output files
				e.write('Time step | Potential Energy | Kinetic Energy| Total Energy|Temperature|Average Temperature\n')
				for i in range(nsteps):
					for j in range(0,Init.ip): #Updating the positions
						vt[j][0]=vt[j][0]+dth*at[j][0]
						vt[j][1]=vt[j][1]+dth*at[j][1]
						vt[j][2]=vt[j][2]+dth*at[j][2]
						rt[j][0]=rt[j][0]+dt*vt[j][0]
						rt[j][1]=rt[j][1]+dt*vt[j][1]
						rt[j][2]=rt[j][2]+dt*vt[j][2]
							
					x=Acceleration.Calc_a(eps,sigma,mass,rt) #Updating the acceleration
					at=x
					for n in range(0, Init.ip): # Updating the velocity
						vt[n]=[vt[n][o]+at[n][o]*dth for o in range(0,3)]

					E_potential[i]=LJ.Calc_LJ(eps,sigma,rt)[3]
					E_kinetic[i]=Kinetic_Energy.Calc_KE(eps,sigma,mass,vt)	
					E_total[i]=E_potential[i]+E_kinetic[i]
					T[i]=(E_kinetic[i]*2)/(3*Init.ip)
					T_sum+=T[i]
					T_avg=T_sum/(i+1) #Calculating running average for temperature
					print ('%f %f %f %f %f %f\n' % (i*dt, E_potential[i], E_kinetic[i], E_total[i],T[i], T_avg))
					e.write(' %f %f %f %f %f%f\n' % (i*dt,E_potential[i],E_kinetic[i],E_total[i], T[i], T_avg))
					f.write('Time step: %f, rx ry rz vx vy vz \n' % (i*dt))
					d.write("%d\nInitial Conditions\n" % (Init.N))					
					for p in range (0, Init.ip):
						f.write('{:d},{:<20.6f},{:<20.6f},{:<20.6f},{:<20.6f},{:<20.6f},{:<20.6f},{:<20.6f},{:<20.6f},{:<20.6f}\n'.format(p, rt[p][0],rt[p][1],rt[p][2],vt[p][0],vt[p][1],vt[p][2],at[p][0],at[p][1],at[p][2]))
						d.write('%s%8.3f%8.3f%8.3f\n' % (Init.Atom_Type,rt[p][0], rt[p][1], rt[p][2]))
				d.close()		
				e.close()
				f.close()
	pylab.plot(E_potential,label='Potential Energy')
	pylab.plot(E_kinetic,label='Kinetic Energy')
	pylab.plot(E_total, label='Total Energy')
	pylab.xlabel('Time')
	pylab.ylabel('Energy')
	pylab.title('Energy Plots, N=%d, L=%f' % (Init.ip, Init.L))
	pylab.savefig('N=%dL=%f_LeapFrog.png' % (Init.ip, Init.L))
	pylab.legend()
	pylab.show()

#Uncomment the next line if you want to specifically run this file
#Run_Leap_Frog(1.0,1.0,1.0,1.0,0.001)
