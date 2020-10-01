#Salman Bin Kashif
#Sarupria Group
#Clemson University
#March 15, 2018

import numpy as np
import Init
import Acceleration
import LJ
import Kinetic_Energy
import pylab
import math

def Run_Vel_ver(eps,sigma,mass,time,dt):
	global rt
	global vt
	global at
	global E_potential
	global E_kinetic
	global T_avg

	Init.initialize()
	Init.print_init()
	

	rt=Init.r #Assigning initial values to poisitons
	vt=Init.v #Assigning initial values to velocities
	at=[[0.0]*3 for i in range(0,Init.ip)] #Initializing the acceleration varible
	sum_a=[[0.0]*3 for i in range(0,Init.ip)] #Initializing the variable to store the a(t) and a(t+dt) values
	nsteps=(int)(time/dt) #Number of steps for simulation to run

	#Initializing the energy and temperature variables

	E_potential=[0.0]*(nsteps)
	E_kinetic=[0.0]*(nsteps)
	E_total=[0.0]*(nsteps)

	at=Acceleration.Calc_a(eps,sigma,mass,rt) #Assigning initial values to acceleration
	
	E_kin=[0.0]*(nsteps)
	T=[0.0]*(nsteps)
	T_avg=0.0
	T_sum=0.0
	print('Time step | Potential Energy | Kinetic Energy| Total Energy| Temperature|Average Temperture\n')

	with open('Vel_verlet_NVE.out','w') as f:
		with open('Vel_verlet_Energy_NVE.out','w') as e:
			e.write('Time step | Potential Energy | Kinetic Energy| Total Energy| Temperature|Average Temperture\n')
			with open('VMD_pos_Vel_ver_NVE.xyz','w') as d:
				for i in range(nsteps): # Loop running in time steps
					for j in range(0,Init.ip): #Updating the positions of all the particles
						rt[j][0]=rt[j][0]+vt[j][0]*dt+0.5*(dt**2)*at[j][0]
						rt[j][1]=rt[j][1]+vt[j][1]*dt+0.5*(dt**2)*at[j][1]
						rt[j][2]=rt[j][2]+vt[j][2]*dt+0.5*(dt**2)*at[j][2]
						
					x=Acceleration.Calc_a(eps,sigma,mass,rt) #Updating the acceleration
					
		
					for l in range(0, Init.ip): #Calculating the sum of a(t) and a(t+delta(t)) to be used in velocity calculation
						sum_a[l][0]=at[l][0]+x[l][0]
						sum_a[l][1]=at[l][1]+x[l][1]
						sum_a[l][2]=at[l][2]+x[l][2]
					at=x
					for n in range(0, Init.ip): # Updating the velocity
						vt[n][0]=vt[n][0]+0.5* dt * (sum_a[n][0])
						vt[n][1]=vt[n][1]+0.5* dt * (sum_a[n][1])
						vt[n][2]=vt[n][2]+0.5*dt * (sum_a[n][2]) 
					E_potential[i]=LJ.Calc_LJ(eps,sigma,rt)[3]
					E_kinetic[i]=Kinetic_Energy.Calc_KE(eps,sigma,mass,vt)	
					E_total[i]=E_potential[i]+E_kinetic[i]
					T[i]=(E_kinetic[i]*2)/(3*Init.ip)
					T_sum+=T[i]
					T_avg=T_sum/(i+1) #Calculating the running average of the temperature
					print('%f %f %f %f %f %f\n' % (i*dt, E_potential[i], E_kinetic[i], E_total[i], T[i], T_avg))
					e.write('%f %20.6f %20.6f %20.6f %20.6f %20.6f\n' % (i*dt,E_potential[i],E_kinetic[i],E_total[i], T[i], T_avg))
					f.write('Time step: %f, rx ry rz vx vy vz \n' % (i*dt))
					d.write("%d\nInitial Conditions\n" % (Init.N))					
					for p in range (0, Init.ip):
						f.write('{:<20.6f},{:<20.6f},{:<20.6f},{:<20.6f},{:<20.6f},{:<20.6f},{:<20.6f},{:<20.6f},{:<20.6f}\n'.format(rt[p][0],rt[p][1],rt[p][2],vt[p][0],vt[p][1],vt[p][2],at[p][0],at[p][1],at[p][2]))
						d.write('%s%8.3f%8.3f%8.3f\n' % (Init.Atom_Type,rt[p][0], rt[p][1], rt[p][2]))
			d.close()		
		e.close()
	f.close()
	pylab.plot(E_potential,label='Potential Energy')
	pylab.plot(E_kinetic,label='Kinetic Energy')
	pylab.plot(E_total, label='Total Energy')
	pylab.xlabel('Time')
	pylab.ylabel('Energy')
	pylab.title('Energy Plots, N=%d, L=%f' % ((Init.ip), Init.L))
	pylab.savefig('N=%dL=%f_Vel_verlet.png' % ((Init.ip), Init.L))
	pylab.legend()
	pylab.show()
			
#Run_Vel_ver(1,1,1,1,0.001)
