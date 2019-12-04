#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --partition=standard
./../../build/main 10 38 80000 100 80000 1 0.5 3.8 700