#This file uses the directory created by gen_jobs and plots SCDWvsT and autocorrelations


import sys
import os
import numpy as np

# Parameters~~~~~~~~~~~~~~~~~~~~~
g = 1    #must be single value
w = 0.5  #must be single value
# Parameters~~~~~~~~~~~~~~~~~~~~~

here = os.path.dirname(os.path.realpath(__file__))
dire = 'g'+str(g)+'w'+str(w)

if not os.path.exists(os.path.join(here,dire)):
	print("The directory corresponding to the input g and w does not exist")
	sys.exit()


subfolders = [f.path for f in os.scandir(os.path.join(here,dire)) if f.is_dir() ]


for folder in subfolders:
	for file in folder:
		if file.endswith(".txt"):
			print(os.path.join("/mydir", file))


	subdir = 'N'+str(N)
	
	if not os.path.exists(os.path.join(here,dire,subdir)):
		os.makedirs(os.path.join(here,dire,subdir))
	for L in La:
		f2 = open(os.path.join(here,dire,subdir,'L'+str(L)+'.sh'),'w+')
		f2.write('#!/bin/bash\n')
		f2.write('#SBATCH --nodes=1\n')
		f2.write('#SBATCH --time=48:00:00\n')
		f2.write('#SBATCH --partition=standard\n')
		f2.write('./../../build/main '+str(N)+' '+str(L)+' 80000 '+str(N*N)+' 80000 '+str(g)+' '+str(w)+' '+str(np.round(L/10.0,1))+' 700')
		
		f3.write('cd N'+str(N)+'\n')
		f3.write('sbatch L'+str(L)+'.sh\n')
		f3.write('cd ..\n')
