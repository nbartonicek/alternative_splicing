#!/bin/bash

module load trimgalore/0.4.4
module load fastqc/0.11.5

############## directory hierarchy ##############
#raw files directory
homedir="/researchers/nenad.bartonicek"

project="PRMT1"
inType="trimgalore"

project_dir="$homedir/projects/$project"
raw_dir="$project_dir/raw"
script_dir="$homedir/projects/$project/scripts"
out_dir="$project_dir/results/"$inType

#extension of the files to be used
inExt="fastq.gz"

#log and command files for bsub
log_dir=$script_dir/"logs"

#make the directory structure   
mkdir -p $out_dir
mkdir -p $log_dir


files=( $(ls $raw_dir/*R1_001.fastq.gz) )
for inFile1 in ${files[@]};do

        inFile2=`echo $inFile1 | sed s/_R1_001.fastq.gz/_R2_001.fastq.gz/`
        
        uniqueID=`basename $inFile1 | sed s/_R1_001.fastq.gz//`
        echo $uniqueID
        
        sample_out_dir=$out_dir/$uniqueID/
        mkdir -p $sample_out_dir
        #echo $command_line

        command_line="trim_galore \
                $inFile1 $inFile2 \
                --gzip \
                --fastqc \
                --paired \
                --quality 20 \
                --length 16 \
                -o $sample_out_dir"
        #echo $command_line
        sbatch -p prod -J trimgalore -t 0-3:00 -c 4 --mem 12G -o $log_dir/%j.out -e $log_dir/%j.err --mail-user=nenad.bartonicek@petermac.org --wrap="$command_line"

done;


