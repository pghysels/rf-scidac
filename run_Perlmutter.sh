#!/bin/bash
#SBATCH -A m3894_g
#SBATCH -C gpu
#SBATCH -q regular
#SBATCH -t 0:30:00
#SBATCH -N 32
#SBATCH -n 128
#SBATCH --ntasks-per-node=4
#SBATCH -c 32
#SBATCH --gpus-per-task=1
#SBATCH --gpu-bind=map_gpu:0,1,2,3


export SLURM_CPU_BIND="cores"
export OMP_NUM_THREADS=16
export MPICH_GPU_SUPPORT_ENABLED=1


srun ./build/main /pscratch/sd/p/pghysels/nstxu_180degree_antenna_phasing_write_matrix_order3/matrix. --sp_compression NONE --sp_enable_gpu > out/order3_Perlmutter_N32_GPU_NONE.log
# srun ./build/main /pscratch/sd/p/pghysels/nstxu_180degree_antenna_phasing_write_matrix_order3/matrix. --sp_compression BLR --sp_disable_gpu > out/order3_Perlmutter_N32_BLR.log


# srun ./build/main /pscratch/sd/p/pghysels/nstxu_180degree_antenna_phasing_write_matrix_order3/matrix. --sp_compression LOSSY --sp_disable_gpu > out/order3_Perlmutter_N32_LOSSY.log
