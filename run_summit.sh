#!/bin/bash
# Begin LSF Directives
#BSUB -P CSC289
#BSUB -W 00:30
#BSUB -nnodes 64
#BSUB -J rf-scidac
#BSUB -o rf-scidac.%J
#BSUB -e rf-scidac.%J

module unload darshan-runtime

export OMP_NUM_THREADS=7

# 64 nodes x 42 cores = 2688 cores total
# every node has 6 GPU's, use 6 resource sets per node (-r 6)
# 1 GPU per resource set (-g 1)
# use 7 threads (cores) per resource set (-c 7  -d packed -b packed:7)
# 7 cores x 6 resource sets = 42 cores per node
# 6 resource sets per node x 64 nodes = 384 MPI ranks
jsrun -n 384 -a 1 -c 7 -g 1 -r 6 -l CPU-CPU -d packed -b packed:7 \
      ./build/main \
      /gpfs/alpine/scratch/pghysels/csc289/rf-scidac/matrix. \
      /gpfs/alpine/scratch/pghysels/csc289/rf-scidac/rhs_0.000000 \
      --sp_compression none \
      --sp_matching 0 \
      --sp_enable_gpu \
      --sp_reordering_method parmetis \
      --help > out/order3_summit_N64_none_parmetis_GPU.log


# --blr_rel_tol 1e-4 \
# --hodlr_butterfly_levels 100 \
# --hodlr_rel_tol 1e-4 \
# --sp_hodlr_min_sep_size 10000000 \
# --sp_print_compressed_front_stats \

# --blr_factor_algorithm COLWISE \
