#!/bin/bash

#module load gi/star/2.3.0e
module load samtools/1.17
module load star/2.5.3a

numcores=15

############## directory hierarchy ##############
#raw files directory
homedir="/researchers/nenad.bartonicek"

project="PRMT1"
in_type="trimgalore"
out_type="star"

project_dir="$homedir/projects/$project"
raw_dir="$project_dir/raw"
script_dir="$homedir/projects/$project/scripts"
in_dir="$project_dir/results/"$in_type
out_dir="$project_dir/results/"$out_type

#extension of the files to be used
inExt="fastq.gz"

#log and command files for bsub
log_dir=$script_dir/"logs/star"

#make the directory structure   
mkdir -p $out_dir
mkdir -p $log_dir


genomeDir=$project_dir"/annotation/star"
#genomeDir="/share/ClusterShare/biodata/contrib/nenbar/genomes/star/mm10_sequin"
#extension of the files to be used
inExt="fq.gz"


samples=`ls $in_dir | sed 's/\\///'` 

for sample in ${samples[@]}; do

        fastq_dir=$in_dir/$sample
	star_dir=$out_dir/$sample

        inFile1=`ls $fastq_dir | grep .*_val_1.$inExt`      
	inFile2=`echo $inFile1 | sed 's/_R1_001_val_1.fq.gz/_R2_001_val_2.fq.gz/'`
        
        star_dir=$out_dir/$sample/
	mkdir -p $star_dir


	starJobName="star."$sample
	indexJobName="index."$sample
	outSortedBam=$outDir"$sample.Aligned.sortedByCoord.out.bam"
	
	star_line="STAR \
                --runMode alignReads \
		--genomeDir $genomeDir \
		--readFilesIn $in_dir/$sample/$inFile1 $in_dir/$sample/$inFile2 \
                --outFileNamePrefix $out_dir/$sample \
		--runThreadN $numcores \
                --outBAMsortingThreadN $numcores \
		--readFilesCommand zcat \
                --twopassMode Basic \
                --outFilterType BySJout \
                --outSAMattributes NH HI AS NM MD\
                --outFilterMultimapNmax 20 \
                --outFilterMismatchNmax 999 \
                --outFilterMismatchNoverReadLmax 0.04 \
                --twopass1readsN -1 \
                --alignEndsType Local \
                --alignSJDBoverhangMin 8 \
                --outSAMtype BAM SortedByCoordinate \
                --outFilterMatchNmin 40"


	index_line="samtools index $outSortedBam"

        sbatch -p prod -J star -t 0-3:00 -c 15 --mem 90G -o $log_dir/%j.out -e $log_dir/%j.err --mail-user=nenad.bartonicek@petermac.org --wrap="$star_line" | cut -d " " -f 4

        #sbatch --dependency=afterok:$job1id -p prod -J star -t 0-2:00 -c 15 --mem 16G -o $log_dir/%j.out -e $log_dir/%j.err --mail-user=nenad.bartonicek@petermac.org --wrap="$index_line"

done;
