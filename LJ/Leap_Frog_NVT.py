#Salman Bin Kashif
#Sarupria Group
#Clemson University
#March 20, 2018

import numpy as np
import Init
import Acceleration
import LJ
import Kinetic_Energy
import matplotlib.pylab as plt
import math

#Leap Frog algorithm for the NVT thermostat

def Run_Leap_Frog(eps,sigma,mass,time,dt,Tn,tau):
	global rt
	global vt
	global at
	global E_potential
	global E_kinetic

	Init.initialize()
	Init.print_init()

	rt=Init.r #Assigning inital values to positions
	vt=Init.v #Assigning initial values to velocity
	at=[[0.0]*3 for i in range(0,Init.ip)]
	sum_a=[[0.0]*3 for i in range(0,Init.ip)]
	dth=dt*0.5
	nsteps=(int)(time/dt)
	
	E_potential=[0.0]*nsteps
	E_kinetic=[0.0]*nsteps
	E_total=[0.0]*nsteps
	E_sum=0.0
	E_avg=0.0

	at=Acceleration.Calc_a(eps,sigma,mass,rt) #Assigning initial values to accelration
	E_kin=[0.0]*(nsteps)
	T=[0.0]*(nsteps)
	
	print('Time step | Potential Energy | Kinetic Energy| Total Energy|Average Total Energy | Temperature\n')
				
	at=Acceleration.Calc_a(eps,sigma,mass,rt) #Assigning initial values to acceleration
	
	with open('Leap_Frog.out','w') as f, open('Leap_Frog_Energy.out','w') as e, open('VMD_pos_leapfrog.xyz','w') as d:
				e.write('Time step | Potential Energy | Kinetic Energy| Total Energy|Average Total Energy | Temperature\n')
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
					#Applying the Brendson Thermostat
					E_kin[i]=Kinetic_Energy.Calc_KE(eps,sigma,mass,vt)
					T[i]=(E_kin[i]*2)/(3*Init.ip)
					lam=math.sqrt(1+(dt/tau)*((Tn/T[i])-1))
					#Applying velocity rescaling
					for s in range(0,Init.ip):
						vt[s][0]=vt[s][0]*lam
						vt[s][1]=vt[s][1]*lam
						vt[s][2]=vt[s][2]*lam
					#Updating the energies
					E_potential[i]=LJ.Calc_LJ(eps,sigma,rt)[3]
					E_kinetic[i]=Kinetic_Energy.Calc_KE(eps,sigma,mass,vt)	
					E_total[i]=E_potential[i]+E_kinetic[i]
					E_sum+=E_total[i]
					E_avg=E_sum/(i+1) #Calculating the running average for energy
					print ('%f %f %f %f %f %f\n' % (i*dt, E_potential[i], E_kinetic[i], E_total[i], E_avg,T[i]))
					e.write(' %f %f %f %f %f %f\n' % (i*dt,E_potential[i],E_kinetic[i],E_total[i],E_avg,T[i]))
					f.write('Time step: %f, rx ry rz vx vy vz \n' % (i*dt))
					d.write("%d\nInitial Conditions\n" % (Init.N))					
					for p in range (0, Init.ip):
						f.write('{:d},{:<20.6f},{:<20.6f},{:<20.6f},{:<20.6f},{:<20.6f},{:<20.6f},{:<20.6f},{:<20.6f},{:<20.6f}\n'.format(p, rt[p][0],rt[p][1],rt[p][2],vt[p][0],vt[p][1],vt[p][2],at[p][0],at[p][1],at[p][2]))
						d.write('%s%8.3f%8.3f%8.3f\n' % (Init.Atom_Type,rt[p][0], rt[p][1], rt[p][2]))
				d.close()		
				e.close()
				f.close()
	plt.plot(E_potential,label='Potential Energy')
	plt.plot(E_kinetic,label='Kinetic Energy')
	plt.plot(E_total, label='Total Energy')
	plt.xlabel('Time')
	plt.ylabel('Energy')
	plt.title('Energy plot, N=%d, L=%f' % (Init.ip, Init.L))
	plt.savefig('N=%dL=%f_Energy_LF_NVT.png' % (Init.ip, Init.L))
	plt.show()
	
	plt.plot(T, label='Temperature')
	plt.xlabel('Time')
	plt.ylabel('Temperature')
	plt.title('Energy plot, N=%d, L=%f' % (Init.ip, Init.L))
	plt.savefig('N=%dL=%f_Temperature_LF_NVT.png' % (Init.ip, Init.L))
	plt.legend()
	plt.show()

#Uncomment the next line if you wish to specifically run this file
#Run_Leap_Frog(1,1,1,1,0.001,0.838038,0.001)
