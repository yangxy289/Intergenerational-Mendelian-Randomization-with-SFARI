#!/bin/bash

input_path="/ASD/SPARK/WGS/5.mr/output/GeneticScore"
output_path="/ASD/SPARK/WGS/5.mr/output/ADD"
genotype_path="/ASD/SPARK/SPARK_phasing_fam/SPARK_WGS_phasing_vcfsample_trio5712_IID_gender_asd.txt"
pc_path="/ASD/SPARK/WGS/4.PCA/output/SPARK_pcair_5712trio_PC5.txt"
phenotype="trait"

prs_path=$input_path/$phenotype
add_path=$output_path/$phenotype
mkdir -p $add_path

for gs in $prs_path/*.txt; do
	prefix=$(basename $gs .txt)
	echo "processing gs: $gs"

	python intergenerational_mr_trio_gwas_multiPGS.py ADD \
		-g $gs \
		-p $genotype_path > $add_path/${prefix}_genotype_trio.tmp
	echo "done ADD genotype"

	python intergenerational_mr_trio_gwas_multiPGS.py ADD \
                -g $add_path/${prefix}_genotype_trio.tmp \
                -p $pc_path > $add_path/${prefix}_ADD_trio.txt
        echo "done ADD pc"

	rm $add_path/${prefix}_genotype_trio.tmp

done
