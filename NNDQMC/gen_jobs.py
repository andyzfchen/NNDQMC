#This file creates a directory in the location it is ran.
#The directory has a name of the g and w values you input
#This Directory has a subdirectory for each size in the Na array
#As well as a script named submit_scripts.sh
#Each size directory has a script file for each imaginary time length
#submit_scripts goes through each directory and submits the jobs within
#must make submit_scripts executable using chmod -u+x submit_scripts.sh



import os
import numpy as np

# Parameters~~~~~~~~~~~~~~~~~~~~~
g = 1    #must be single value
w = 0.5  #must be single value
Na = np.array([6,8,10,12,14]) #array of lattice sizes
La = np.array([37,38,39,40,41,42,43,44,45,46]) #array of imaginary time slices
# Parameters~~~~~~~~~~~~~~~~~~~~~

here = os.path.dirname(os.path.realpath(__file__))
dire = 'g'+str(g)+'w'+str(w)

if not os.path.exists(os.path.join(here,dire)):
	os.makedirs(os.path.join(here,dire))

f3 = open(os.path.join(here,dire,'submit_scripts.sh'),'w+')


for N in Na:
	subdir = 'N'+str(N)
	if not os.path.exists(os.path.join(here,dire,subdir)):
		os.makedirs(os.path.join(here,dire,subdir))
	for L in La:
		f2 = open(os.path.join(here,dire,subdir,'L'+str(L)+'.sh'),'w+')
		f2.write('#!/bin/bash\n')
		f2.write('#SBATCH --nodes=1\n')
		if(N==6):
			f2.write('#SBATCH --time=02:00:00\n')
		if(N==8):
			f2.write('#SBATCH --time=04:00:00\n')
		if(N==10):
			f2.write('#SBATCH --time=10:00:00\n')
		if(N==12):
			f2.write('#SBATCH --time=20:00:00\n')
		if(N==14):
			f2.write('#SBATCH --time=30:00:00\n')
		f2.write('#SBATCH --partition=standard\n')
		f2.write('./../../build/main '+str(N)+' '+str(L)+' 80000 '+str(N*N)+' 80000 '+str(g)+' '+str(w)+' '+str(np.round(L/10.0,1))+' 700')
		
		f3.write('cd N'+str(N)+'\n')
		f3.write('sbatch L'+str(L)+'.sh\n')
		f3.write('cd ..\n')
