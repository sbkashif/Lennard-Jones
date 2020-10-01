#Salman Bin Kashif
#Sarupria Group
#Clemson University
#Date created:
#Last modified:

import numpy as np
import Init
import Acceleration
import LJ
import Kinetic_Energy
import pylab
import math

def Run_DPD(eps,sigma,mass,time,dt,Tn):
	global rt
	global vt
	global at
	global E_potential
	global E_kinetic

#Initialization
	Init.initialize()
	Init.print_init()

	rt=Init.r
	vt=Init.v

	at=[[0.0]*3 for i in range(0,Init.ip)]
	nsteps=(int)(time/dt)
	
	E_potential=[0.0]*(nsteps)
	E_kinetic=[0.0]*(nsteps)
	E_total=[0.0]*(nsteps)
	T=[0.0]*nsteps

	gamma=1.0
	sig=np.sqrt(2*gamma*Tn)

	at=Acceleration.Calc_a(eps,sigma,mass,rt,vt,sig,gamma,dt)
	
	E_sum=0.0
	E_avg=0.0

	print('Time step | Potential Energy | Kinetic Energy | Total Energy | Average Total Energy | Temperature\n')

	with open('DPD_NVT.out','w') as f:
		with open('DPD_Energy_NVT.out','w') as e:
			e.write('Time step | Potential Energy | Kinetic Energy | Total Energy | Average Total Energy | Temperature\n')
			with open('VMD_pos_NVT.xyz','w') as d:
				for i in range(nsteps):
					for j in range(0, Init.ip):
						vt[j][0]=vt[j][0]+0.5*dt*at[j][0]
						vt[j][1]=vt[j][1]+0.5*dt*at[j][1]
						vt[j][2]=vt[j][2]+0.5*dt*at[j][2]
						rt[j][0]=rt[j][0]+dt*vt[j][0]
						rt[j][1]=rt[j][1]+dt*vt[j][1]
						rt[j][2]=rt[j][2]+dt*vt[j][2]

					at=Acceleration.Calc_a(eps,sigma,mass,rt,vt,sig,gamma,dt)
					for s in range(0, Init.ip):
						vt[s][0]=vt[s][0]+0.5*dt*at[s][0]
						vt[s][1]=vt[s][1]+0.5*dt*at[s][1]
						vt[s][2]=vt[s][2]+0.5*dt*at[s][2]
					
					#at=Acceleration.Calc_a(eps,sigma,mass,rt,vt,sig,gamma,dt)
					
					E_potential[i]=LJ.Calc_LJ(eps,sigma,rt,vt,sig,gamma,dt)[3]
					E_kinetic[i]=Kinetic_Energy.Calc_KE(eps,sigma,mass,vt)
					T[i]=(E_kinetic[i]*2)/(3*Init.ip)
					E_total[i]=E_potential[i]+E_kinetic[i]
					E_sum+=E_total[i]
					E_avg=E_sum/(i+1)				
					print('%f %f %f %f %f %f\n' % (i*dt, E_potential[i], E_kinetic[i], E_total[i],E_avg,T[i]))
                                        e.write('%f %f %f %f %f %f\n' % (i*dt,E_potential[i],E_kinetic[i],E_total[i],E_avg,T[i]))
                                        f.write('Time step: %f, rx ry rz vx vy vz \n' % (i*dt))
                                        d.write("%d\nInitial Conditions\n" % (Init.N))
                                        for p in range (0, Init.ip):
	                                        f.write('{:<20.6f},{:<20.6f},{:<20.6f},{:<20.6f},{:<20.6f},{:<20.6f},{:<20.6f},{:<20.6f},{:<20.6f}\n'.format(rt[p][0],rt[p][1],    rt[p][2],vt[p][0],vt[p][1],vt[p][2],at[p][0],at[p][1],at[p][2]))
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
        pylab.savefig('N=%dL=%f_DPD_NVT.png' % (Init.ip, Init.L))

        pylab.plot(T, label='Temperature')
        pylab.xlabel('Time')
        pylab.ylabel('Temperature')
        pylab.title('Temperature plot, N-%d, L=%f' % (Init.ip, Init.L))
        pylab.legend()
        pylab.show()
        pylab.savefig('N=%dL=%f_DPD_NVT.png' % (Init.ip, Init.L))
 
#Uncomment the next line if you want to specifically run this file                      
Run_DPD(1,1,1,10,0.001,0.838038)
