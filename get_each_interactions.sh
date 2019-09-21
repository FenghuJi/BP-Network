#!/usr/bin/bash
work_dir=/home/jfh/data/Virome/Gajewski
filelist=$work_dir/filelist
mkdir -p $work_dir/08_interactions_all
cat $filelist|while read LINE;
do
	sample=$LINE
	cat $work_dir/07_interaction/prophage_all_interaction.txt|grep $sample|grep species|cut -f1,4 >  $work_dir/08_interactions_all/$sample'.csv'
	cat $work_dir/06_crispr/06_interaction/crispr_interactions.txt|grep $sample|grep species|cut -f1,4 >>  $work_dir/08_interactions_all/$sample'.csv'
done
sed -i 's/\t/,/' $work_dir/08_interactions_all/*
