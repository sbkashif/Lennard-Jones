***Molecular Dynamics cum Monte Carlo Simulator***

#Developed by: Salman Bin Kashif
#Affiliation: Sarupria Group, Clemson University, South Carolina
#Date: March 23, 2018

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

1) Acceleration.py
2) Init.py
3) Kinetic_Energy.py
4) Leap_Frog_NVE.py
5) Leap_Frog_NVT.py
6) LJ_MC.py
7) LJ.py
8) MC.py
9) reset.sh
10) Run.py
11) Vel_verlet_NVE.py
12) Vel_verlet_NVT.py

**Run the code**

#Enter the following command:

python2 Run.py

#Follow the instructions on screen and enter the required input values. Please pay due care to the warnings after the input values have been added. The input values that you have enetered may possibly be against a realistic setup. 

# Once the run has ended:

1) Copy the needed output files to another folder

2) Enter the following command: 

bash reset.sh

* This reset command is basically clearling up the output files of previous run and leaving a clean slate for the new run. The code is designed to replace the old files in new run but it is recommended to reset.

Note:-

After any run or at any point, if the user wish to visualize the system, type the following command:

vmd <filename.xyz>

*This file is generated as a part of the code.