#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=30:00:00
#SBATCH --partition=standard
./../../build/main 14 43 80000 196 80000 1 0.5 4.3 700