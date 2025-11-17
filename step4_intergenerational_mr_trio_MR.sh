#!/bin/bash

input_path="/ASD/SPARK/WGS/5.mr/output/ADD"
output_path="/ASD/SPARK/WGS/5.mr/output/MR"
phenotype="trait"

add_path=$input_path/$phenotype/
mr_path=$output_path/$phenotype/trio5712_intergenerationalmr
mkdir -p $mr_path

for gs in $add_path/*.txt; do
	prefix=$(basename $gs _ADD_trio.txt)
	echo "processing: $gs"

	python intergenerational_mr_trio_gwas_multiPGS.py MR \
		-I $gs \
		-x h1,h2,h3,h4 \
		--covar gender,PC1,PC2,PC3,PC4,PC5 \
		-y pheno > $mr_path/${prefix}_MR_trio_5712_h14_gender_PC.txt
	echo "done h14 MR"

done