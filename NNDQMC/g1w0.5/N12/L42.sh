#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=20:00:00
#SBATCH --partition=standard
./../../build/main 12 42 80000 144 80000 1 0.5 4.2 700