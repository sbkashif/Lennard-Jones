#Salman Bin Kashif
#Sarupria Research Group
#Clemson University
#Started: April 16, 2018
#Last modified:



import numpy as np
import Init
import Acceleration
import LJ
import Kinetic_Energy
import pylab
import math
 

def dW1(delta_t):
	return np.random.normal(loc=0.0,scale=np.sqrt(delta_t))

def dW3(delta_t):
	return np.random.normal(loc=0.0,scale=np.sqrt(delta_t**(3.0)/3.0))

def Run_Langevin(eps,sigma,mass,time,dt,Tn):
       global rt
       global vt
       global at
       global E_potential
       global E_kinetic

#Initialization

       Init.initialize()
       Init.print_init()


       rt=Init.r #Assigning initial values to poisitons
       vt=Init.v #Assigning initial values to velocities

       at=[[0.0]*3 for i in range(0,Init.ip)]
       sum_a=[[0.0]*3 for i in range(0,Init.ip)]
       nsteps=(int)(time/dt)
       

       E_potential=[0.0]*(nsteps)
       E_kinetic=[0.0]*(nsteps)
       E_total=[0.0]*(nsteps)

       at=Acceleration.Calc_a(eps,sigma,mass,rt) #Assigning initial values to accelration
       E_kin=[0.0]*(nsteps)
       T=[0.0]*(nsteps)
       E_sum=0.0
       E_avg=0.0
       gamma=1.0
       #E_k=Kinetic_Energy.Calc_KE(eps,sigma,mass,vt)
       #T0=(E_k*2)/(3*Init.ip)
       sig=np.sqrt(2*gamma*Tn)
       print('Time step | Potential Energy | Kinetic Energy| Total Energy| Average Total Energy| Temperature\n')

	
       with open('Vel_verlet_NVT.out','w') as f, open('Vel_verlet_Energy_NVT.out','w') as e:
                       e.write('Time step | Potential Energy | Kinetic Energy| Total Energy| Average Total Energy| Temperature\n')
                       with open('VMD_pos_NVT.xyz','w') as d:
                               	for i in range(nsteps): # Loop running in time steps

					for j in range(0,Init.ip): #Updating the positions of all the particles
                                       		rt[j][0]=rt[j][0]+vt[j][0]*dt+0.5*(dt**2)*(at[j][0]-gamma*vt[j][0])+sig*dW3(dt)
                                       		rt[j][1]=rt[j][1]+vt[j][1]*dt+0.5*(dt**2)*(at[j][1]-gamma*vt[j][1])+sig*dW3(dt)
                                       		rt[j][2]=rt[j][2]+vt[j][2]*dt+0.5*(dt**2)*(at[j][2]-gamma*vt[j][2])+sig*dW3(dt)

#Storing the updated forces for the next step 
                                	atu=Acceleration.Calc_a(eps,sigma,mass,rt) 

#Calculating the sum of a(t) and a(t+delta(t)) to be used in velocity calculation

                                 	for m in range(0, Init.ip):
                                        	sum_a[m][0]=at[m][0]+atu[m][0]
                                        	sum_a[m][1]=at[m][1]+atu[m][1]
                                        	sum_a[m][2]=at[m][2]+atu[m][2]
#Updating the velocity   

	                         	for n in range(0, Init.ip): # Updating the velocity
                                       		vt[n][0]=vt[n][0]+0.5*dt*(sum_a[n][0])-(gamma*vt[n][0]*dt)+sig*dW1(dt)-gamma*(0.5*(dt**2)*(at[n][0]-gamma*vt[n][0])+sig*dW3(dt))
                                       		vt[n][1]=vt[n][1]+0.5*dt*(sum_a[n][1])-(gamma*vt[n][1]*dt)+sig*dW1(dt)-gamma*(0.5*(dt**2)*(at[n][1]-gamma*vt[n][1])+sig*dW3(dt))
						vt[n][2]=vt[n][2]+0.5*dt*(sum_a[n][2])-(gamma*vt[n][2]*dt)+sig*dW1(dt)-gamma*(0.5*(dt**2)*(at[n][2]-gamma*vt[n][2])+sig*dW3(dt))
                                 	at=atu
#Updating the energies
                                 	E_potential[i]=LJ.Calc_LJ(eps,sigma,rt)[3]
                                 	E_kinetic[i]=Kinetic_Energy.Calc_KE(eps,sigma,mass,vt)
                                 	E_total[i]=E_potential[i]+E_kinetic[i]
					T[i]=(E_kinetic[i]*2)/(3*Init.ip)
                                 	E_sum+=E_total[i]
                                 	E_avg=E_sum/(i+1)
                                 	print('%f %f %f %f %f %f\n' % (i*dt, E_potential[i], E_kinetic[i], E_total[i],E_avg,T[i]))
                                       	e.write('%f %f %f %f %f %f\n' % (i*dt,E_potential[i],E_kinetic[i],E_total[i],E_avg,T[i]))
                                       	f.write('Time step: %f, rx ry rz vx vy vz \n' % (i*dt))
                                       	d.write("%d\nInitial Conditions\n" % (Init.N))
                                       	for p in range (0, Init.ip):
                                               f.write('{:<20.6f}\t{:<20.6f}\t{:<20.6f}\t{:<20.6f}\t{:<20.6f}\t{:<20.6f}\t{:<20.6f}\t{:<20.6f}\t{:<20.6f}\n'.format(rt[p][0],rt[p][1],rt[p][2],vt[p][0],vt[p][1]    ,vt[p][2],at[p][0],at[p][1],at[p][2]))
                                               d.write('%s%8.3f%8.3f%8.3f\n' % (Init.Atom_Type,rt[p][0], rt[p][1], rt[p][2]))
                       d.close()
       e.close()
       f.close()
       pylab.plot(E_potential,label='Potential Energy')
       pylab.plot(E_kinetic,label='Kinetic Energy')
       pylab.plot(E_total, label='Total Energy')
       pylab.xlabel('Time')
       pylab.ylabel('Energy')
       pylab.title('Energy plots, N=%d, L=%f' % (Init.ip, Init.L))
       pylab.legend()
       pylab.show()
       pylab.savefig('N=%dL=%f_VerVer_NVT.png' % (Init.ip, Init.L))
 
       pylab.plot(T, label='Temperature')
       pylab.xlabel('Time')
       pylab.ylabel('Temperature')
       pylab.title('Temperature plot, N-%d, L=%f' % (Init.ip, Init.L))
       pylab.legend()
       pylab.show()
       pylab.savefig('N=%dL=%f_VelVer_NVT.png' % (Init.ip, Init.L))
	
Run_Langevin(1,1,1,10,0.001,0.838038)
