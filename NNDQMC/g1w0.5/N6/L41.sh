#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=02:00:00
#SBATCH --partition=standard
./../../build/main 6 41 80000 36 80000 1 0.5 4.1 700