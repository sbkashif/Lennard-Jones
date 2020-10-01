#Salman Bin Kashif
#Sarupria Group
#Clemson University
#Date started: April 30, 2018
#Last modified: May 2,2018

import numpy as np
import math


def extract_vel(input_file):
	global v
	global v0
	global v0_1
	global start
	global vacf
	global v_avg
	global dth
	global nsteps

	t=10		#Total time
	dt=0.001	#Time step size

	nsteps=(int)(t/dt)
	init=open('Init_Cond.xyz','r')
	N=(int)(init.readlines()[0])
	init.close()
	end=N
	v=[0.0]*3
	v0=[0.0]*3
	vi=[0.0]*3
	init.close()
	init=open('Init_Vel.out','r')
	plc_hold=init.readlines()
	y=plc_hold[0]

	a=[z.strip() for z in y.split (' ')]
	v0[0]=(float)(a[1])
	v0[1]=(float)(a[3])
	v0[2]=(float)(a[4])
	vi[0]=(float)(a[1])
	vi[1]=(float)(a[3])
	vi[2]=(float)(a[4])

	v_avg=[]


	vacf=[]
	print vi



	start=0
	dt=0.001
	t=10.000
	nsteps=(int)(t/dt)
	dth=0.000
	with open(input_file, 'r') as f:
		with open('VACF.out','w') as vc:
			x=f.readlines()
			v_den=np.dot(v0,v0)
			v_den_avg=np.mean(v_den)
			v_num=np.dot(v0,v0)
			v_num_avg=np.mean(v_num)
			v_avg.append(v_num_avg/v_den_avg)
			vacf.append(np.mean(v_avg[0]))
			vc.write('%20.6f\t%20.6f\n'%(0.0000,vacf[0]))

			s=0.001
			i=0.001		
			start_num=1
			start_den=1
			while s<(t-0.001):
				while i < (t-0.001):
					print ('num=%f'%(start_num))
					y=x[start_num]
					a=[z.strip() for z in y.split (',')]
					print('%f\t%d\t%f\n'%(i,start_num,s))
					v[0]=(float)(a[3])
					v[1]=(float)(a[4])
					v[2]=(float)(a[5])		
					v_num=np.dot(v,v0)
					v_den=np.dot(v0,v0)
					y=x[start_den]
					a=[z.strip() for z in y.split(',')]
					v0[0]=(float)(a[3])
					v0[1]=(float)(a[4])
					v0[2]=(float)(a[5])
					v_avg.append(v_num/v_den)
					i=i+s	
					start_den=start_den+N+1
					start_num=(int)(round(start_den+(s/dt)*(N+1)))
					
					
					
					
				vacf.append(np.mean(v_avg))
				print('%20.6f\t%20.6f\n'%((s),vacf[(int)(s*1000)]))
				vc.write('%20.6f\t%20.6f\n'%((s),vacf[(int)(s*1000)]))
				s=s+dt
				i=s
				start_den=1
				print ('s=%f\tdt=%f\tN=%f'%((start_den+(s/dt)*(N+1)),dt,(N+1)))
				print s
				start_num=(start_den+((s/dt)*(N+1)))
				start_num=(int)(round(start_num))
				print ('ls=%d'%(start_num))
				v0=vi
			vc.close()
	f.close()
		
			
extract_vel('DPD_NVT.out')
		
	
