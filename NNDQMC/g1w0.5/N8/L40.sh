#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=04:00:00
#SBATCH --partition=standard
./../../build/main 8 40 80000 64 80000 1 0.5 4.0 700