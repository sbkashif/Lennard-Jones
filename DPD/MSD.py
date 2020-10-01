#Salman Bin Kashif
#Sarupria Group
#Clemson University
#Date started: May 1, 2018
#Last modified: May 2, 2018

import numpy as np

def extract_pos(input_file):

	global v
	global start
	global vacf
	global v_avg
	global dth
	global nsteps

	t=10
	dt=0.001

	nsteps=(int)(t/dt)
	init=open('Init_Cond.xyz','r')
	N=(int)(init.readlines()[0])
	init.close()
	end=N
	r=[0.0]*3
	r0=[0.0]*3
	ri=[0.0]*3
	init=open('Init_Cond.xyz','r')
	plc_hold=init.readlines()[3]
	a=[z.strip() for z in plc_hold.split (' ')]
	print a
	r0[0]=(float)(a[3])
	r0[1]=(float)(a[6])
	r0[2]=(float)(a[9])	
	ri[0]=(float)(a[3])
	ri[1]=(float)(a[6])
	ri[2]=(float)(a[9])

	r_avg=[]


	msd=[]
	
	dt=0.001
	t=10.0000

	with open(input_file, 'r') as f:
		with open('MSD.out','w') as ms:
			x=f.readlines()
			diff=[a_i-b_i for a_i,b_i in zip(r0,r0)]
			diff2_terms = [i*i for i in diff]
			diff2=np.sum(diff2_terms)
			r_avg.append(diff2)
			print r_avg[0]
			msd.append(np.mean(r_avg[0]))
			print msd[0]
			ms.write('%20.6f\t%20.6f\n'%(0.0000,msd[0]))
			
			s=0.001
			i=0.001
			start_term=1
			start_ref=1
			while s< (t-0.001):
				while i< (t-0.001):
					print ('i=%f'%(i))
					y=x[start_term]
					a=[z.strip() for z in y.split (',')]
					r[0]=(float)(a[0])
					r[1]=(float)(a[1])
					r[2]=(float)(a[2])
					diff=[a_i-b_i for a_i,b_i in zip(r,r0)]
					diff2_terms=[term *term for term in diff]
					diff2=np.sum(diff2_terms)
					r_avg.append(diff2)
					y=x[start_ref]
					a=[z.strip() for z in y.split(',')]
					r0[0]=(float)(a[0])
					r0[1]=(float)(a[1])
					r0[2]=(float)(a[2])
					i=i+s	
					start_ref=start_ref+N+1
					start_term=(int)(round(start_ref+(s/dt)*(N+1)))
				msd.append(np.mean(r_avg))
				print('%20.6f\t%20.6f\n'%((s),msd[(int)(s*1000)]))
				ms.write('%20.6f\t%20.6f\n'%((s),msd[(int)(s*1000)]))
				s=s+dt
				i=s
				start_ref=1
				start_term=(start_ref+((s/dt)*(N+1)))
				start_term=(int)(round(start_term))
				ro=ri
			ms.close()
		f.close()
extract_pos('DPD_NVT.out')

