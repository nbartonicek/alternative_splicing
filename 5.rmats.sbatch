#!/bin/bash

#SBATCH -J rmats                          # set name for job
#SBATCH --partition prod                       # set partition (prod, prod_med, prod_short)
#SBATCH --time=0-04:00                         # set walltime 
#SBATCH --cpus-per-task=12                      # no. cpus should match the no. threads used in script e.g. -t below
#SBATCH --mem=96G                               # set total memory for job
#SBATCH -o logs/%j.out                              # set name of stdout file (%j is jobid)
#SBATCH -e logs/%j.err                              # set name of stderr file
#SBATCH --mail-type=ALL                        # notify on job status (ALL is BEGIN, END, FAIL, INVALID_DEPEND, REQUEUE, STAGE_OUT)
#SBATCH --mail-user=nenad.bartonicek@petermac.org   # change to your email address

#module load rmats/4.0.2
 
module load samtools/1.17
module load star/2.7.5b
module load bamtools/2.5.0
module load cmake/3.16.5
module load gsl/2.7.1


############## directory hierarchy ##############
#raw files directory
homedir="/researchers/nenad.bartonicek"

project="PRMT1"
in_type="star"
out_type="rmats"

project_dir="$homedir/projects/$project"
raw_dir="$project_dir/raw"
script_dir="$homedir/projects/$project/scripts"
in_dir="$project_dir/results/"$in_type
out_dir="$project_dir/results/"$out_type
treatment_dir="/researchers/nenad.bartonicek/projects/PRMT1/results/treatments"

#extension of the files to be used
inExt="sortedByCoord.out.bam"

#log and command files for bsub
log_dir=$script_dir/"logs/rmats"

gtfFile="/researchers/nenad.bartonicek/projects/PRMT1/annotation/Mus_musculus.GRCm39.110.gtf"
condition_dir=$out_dir"/rMATS_WTuntreated_Prmt1KOuntreated_take2/"

mkdir -p $condition_dir

nthread=12
stranded=fr-firststrand
readlength=66
condition1=$treatment_dir"/wt_untreated"
condition2=$treatment_dir"/Prmt1ko_untreated"
STARindex="/researchers/nenad.bartonicek/projects/PRMT1/annotation/star"
 
cd $condition_dir
echo $gtfFile
echo $outDir
echo $nthread
echo $stranded
echo $readlength
 
rmats.py --b1 $condition1 \
   --b2 $condition2 \
   --gtf $gtfFile \
   --novelSS \
   --variable-read-length \
   --allow-clipping \
   --bi $STARindex \
   --tmp $condition_dir/tmp \
   --od $condition_dir \
   --libType $stranded \
   --nthread $nthread \
   --readLength $readlength \
   --task both

echo FINISHING
 
