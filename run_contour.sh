#! /bin/bash
#$ -l h_rt=5:00:00
#$ -cwd
#$ -V
#$h_data=2G

module load python/2.7
python c_runner.py