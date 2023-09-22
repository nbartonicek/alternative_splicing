#!/bin/bash

#SBATCH -J rsem_index                          # set name for job
#SBATCH --partition prod                       # set partition (prod, prod_med, prod_short)
#SBATCH --time=0-04:00                         # set walltime 
#SBATCH --cpus-per-task=8                      # no. cpus should match the no. threads used in script e.g. -t below
#SBATCH --mem=32G                               # set total memory for job
#SBATCH -o logs/%j.out                              # set name of stdout file (%j is jobid)
#SBATCH -e logs/%j.err                              # set name of stderr file
#SBATCH --mail-type=ALL                        # notify on job status (ALL is BEGIN, END, FAIL, INVALID_DEPEND, REQUEUE, STAGE_OUT)
#SBATCH --mail-user=nenad.bartonicek@petermac.org   # change to your email address


module purge

module load rsem/1.3.0
module load star/2.5.3a

annotation_dir="/researchers/nenad.bartonicek/projects/PRMT1/annotation"
rsem-prepare-reference --gtf $annotation_dir/Mus_musculus.GRCm39.110.gtf \
               --star \
               $annotation_dir/Mus_musculus.GRCm39.dna.toplevel.fa \
               $annotation_dir/rsem/mm39