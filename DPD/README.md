***Langevin cum DPD simulator with mesocopic functions***

#Developed by: Salman Bin Kashif
#Affiliation: Sarupria Group, Clemson University, South Carolina
#Date: May 2, 2018

#Purpose: Class Assignment: ChE 8450/Multiscale Modelling

***Running Instructions***

**Make sure that Python 2.7 module is loaded in the Plametto server**

#Some of the statements written in this code are valid only for Python 2.7 versions and would not work for Python 3. If there is an error thrown for the print steps or pylab, it is surely coming from Python 2 not loaded.
#Python 2.7 and related libraries can be loaded on the Palmetto server by following steps::

1) module add anaconda3/2.5.0
2) conda create -n my_env python=2.7
3) source activate my_env
4) conda install numpy=1.1
5) conda install matplotlib


**Unzip the folder**

#Please make sure that the following files are present.

1) Langevin
2) DPD

Langevin will have:

1) Acceleration.py
2) Init.py
3) Kinetic_Energy.py
4) LJ.py
5) reset.sh
6) Langevin.py

DPD will have:

1) Acceleration.py
2) Init.py
3) Kinetic_Energy.py
4) LJ.py
5) reset.sh
6) Dpd.py
7) VACF.py
8) MSD.py

**Run the code**

For Langevin:

#Go to the respective folder

#Enter the following command:

python Langevin.py

#Follow the instructions on screen and enter the required input values. Please pay due care to the warnings after the input values have been added. The input values that you have enetered may possibly be against a realistic setup. 

#The code is designed to return the plots. In case it didn't, please go to GNUPLOT

#Type the following command:

p 'Vel_verlet_Energy_NVT.out' u 1:6 w l

#This will give temperature curve

#For DPD:

#Go to the respective folder

#Enter the following command:

python DPD.py

#Follow the instructions on screen and enter the required input values. Please pay due care to the warnings after the input values have been added. The input values that you have enetered may possibly be against a realistic setup. 

#For VACF:

#Make sure you are in DPD folder and the run has ended.
#Make sure the following files are present:

DPD_Energy_NVT.out
DPD_NVT.out
Init_Vel.out
Interact_Force.out
Init_cond.xyz

#Type the following command

python VACF.py

#Once the run has ended.

#Go to GNUPLOT

#Type the following command

p 'VACF.out' u 1:2 w l

This will give the curve

#For MSD:

#Make sure you are in DPD folder and the run has ended.
#Make sure the following files are present:

Vel_verlet_Energy_NVT.out
Vel_verlet_NVT.out
Init_Vel.out
Interact_Force.out
Init_cond.xyz

#Type the following command

python MSD.py

#Once the run has ended.

#Go to GNUPLOT

#Type the following command

p 'MSD.out' u 1:2 w l

# Once the run has ended and if you wish to run again, please reset it:

1) If needed, copy the output files to another folder

2) Enter the following command: 

bash reset.sh

* This reset command is basically clearling up the output files of previous run and leaving a clean slate for the new run. The code is designed to replace the old files in new run but it is recommended to reset.

Note:-

After any run or at any point, if the user wish to visualize the system, type the following command:

vmd <filename.xyz>

*This file is generated as a part of the code.