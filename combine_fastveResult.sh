#！/usr/bin/bash
# bash /home/jfh/data/metagenome/39_fastVE_result/gajewski /home/jfh/data/metagenome/39_file_list
help()
{
    cat <<- EOF
    Desc: The script is used to combine fastVE result and create a table filled with phage abundance
    Usage: bash combine_fastVE.sh <result_dir> <sample_list>
    Author: Hugo Ji
EOF
    exit 0
}

case $1 in
        -h|-help)
                help;;
esac
exit 0

result_dir=$1;
sample_list=$2;
touch 12_phage_abundance.txt
sample1 = =`ls $result_dir |head -n 1`
cat $result_dir/$sample1/abundance.tsv|cut -f1 > 12_phage_abundance.txt

# combine the results
cat $sample_list|while read LINE;
do
	sample = $LINE;
	cat $result_dir/$sample/abundance.tsv|cut -f5 |sed 's/tpm/'$sample'/' > $sample'.abundance.temp'
	paste 12_phage_abundance.txt $sample'.imgvr.abundance.temp'> phage_paste.temp
	mv phage_paste.temp 12_phage_abundance.txt
	rm *.temp
done