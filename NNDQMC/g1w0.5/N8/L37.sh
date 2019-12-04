#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=04:00:00
#SBATCH --partition=standard
./../../build/main 8 37 80000 64 80000 1 0.5 3.7 700