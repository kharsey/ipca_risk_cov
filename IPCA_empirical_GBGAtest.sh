#!/bin/bash
 
#SBATCH -J IPCA_empirical_GBGAtest_K6				# job name
##SBATCH --constraint IvyBridge
#SBATCH -n 4
#SBATCH -p serial
#SBATCH --mem-per-cpu=4000
#SBATCH -t 1-0:00                  # wall time (D-HH:MM)
#SBATCH -o outslurm/slurm.%A.%a.out             # STDOUT (%j = JobId)
#SBATCH -e outslurm/slurm.%A.%a.err             # STDERR (%j = JobId)
#SBATCH --array=1-5

module load matlab/2017a

matlab -nodesktop -nosplash << EOF
addpath /home/sjpruitt/CrossSectionFactors
addpath /home/sjpruitt/Data
addpath /home/sjpruitt/
K = 6;
dataname = 'IPCADATA_FNW36_RNKDMN_CON'
clustersuffix = ['cluster' num2str($SLURM_ARRAY_TASK_ID) 'of5'];
bootsims = 1000;
run('IPCA_empirical_GBGAtest.m')
EOF

