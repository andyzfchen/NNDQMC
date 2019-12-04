#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=04:00:00
#SBATCH --partition=standard
./../../build/main 8 43 80000 64 80000 1 0.5 4.3 700