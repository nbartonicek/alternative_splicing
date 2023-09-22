#!/bin/bash

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
condition_dir=$out_dir"/rMATS_WT_untreated_Prmt1KO_untreated/"

mkdir -p $condition_dir

nthread=12
stranded=fr-unstranded
readlength=66
condition1=$treatment_dir"/wt_untreated"
condition2=$treatment_dir"/Prmt1ko_untreated"
condition3=$treatment_dir"/wt_treated"
condition4=$treatment_dir"/Prmt1ko_treated"

STARindex="/researchers/nenad.bartonicek/projects/PRMT1/annotation/star"
 
cd $condition_dir
echo $gtfFile
echo $outDir
echo $nthread
echo $stranded
echo $readlength

analyses=("$condition1 $condition2" \
"$condition3 $condition4" \
"$condition1 $condition3" \
"$condition2 $condition4")

for analysis in "${analyses[@]}"; do
	set -- $analysis
	#echo "$1 vs $2"
	c1=`basename $1`
	c2=`basename $2`
	condition_dir=$out_dir/$c1"_"$c2/
	echo $condition_dir
	mkdir -p $condition_dir
	sbatch -p prod -J rmats$c1$c2 -t 0-3:00 --cpus-per-task=12 --mem 24G -o $log_dir/%j.out -e $log_dir/%j.err --mail-user=nenad.bartonicek@petermac.org \
	--wrap="rmats.py --b1 $treatment_dir/$c1 \
	   --b2 $treatment_dir/$c2 \
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
	   --task both"
done 
 
echo FINISHING
 
