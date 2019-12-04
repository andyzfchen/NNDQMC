#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=30:00:00
#SBATCH --partition=standard
./../../build/main 14 37 80000 196 80000 1 0.5 3.7 700