# BP-Network
scripts for bacteria &amp; phage network constructione   
assembly_reads.pl : 将不同的phage来源组合在一起并修改名字，用于后续的聚类
extract_flanking_seq.py： 提取prophage两端的flanking序列   
get_cut_crispr.pl：获取crispr阵列两端的flanking序列   
get_spacer.pl ：从crispr预测结果中提取出spacer的序列   
prophage_taxaon_filter.py ：对 propage的flanking结果处理后的比对文本初步过滤，找到每个contig的最优分类结果   
specident_CUT_BLAST.pl ：解析blast结果，获得每个flanking序列的blast，包括LCA的计算，最优比对结果以及coverage过滤等   
prophage_PC2host.py： 获取PC id和bacteria id的互作表格   
