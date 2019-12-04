#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --partition=standard
./../../build/main 10 40 80000 100 80000 1 0.5 4.0 700