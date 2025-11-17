#!/bin/bash

input_path="/ASD/SPARK/WGS/5.mr/output/ADD"
output_path="/ASD/SPARK/WGS/5.mr/output/MR"
phenotype="trait"

add_path=$input_path/$phenotype
mr_path=$output_path/$phenotype/trio5712_gt_regression
mkdir -p $mr_path

# use genotype_asd_regression function

for gs in $add_path/*.txt; do
	prefix=$(basename $gs _ADD_trio.txt)
	echo "processing: $gs"

	python intergenerational_mr_trio_gwas_multiPGS.py MR \
                -I $gs \
                -x child_genotype_score \
                --covar gender,PC1,PC2,PC3,PC4,PC5 \
                -y pheno > $mr_path/${prefix}_trio5712_gt_gender_PC.txt
        echo "done gt MR"

	python intergenerational_mr_trio_gwas_multiPGS.py MR \
                -I $gs \
                -x father_genotype_score \
                --covar gender,PC1,PC2,PC3,PC4,PC5 \
                -y pheno >> $mr_path/${prefix}_trio5712_gt_gender_PC.txt
        echo "done gt MR"

	python intergenerational_mr_trio_gwas_multiPGS.py MR \
                -I $gs \
                -x maternal_genotype_score \
                --covar gender,PC1,PC2,PC3,PC4,PC5 \
                -y pheno >> $mr_path/${prefix}_trio5712_gt_gender_PC.txt
        echo "done gt MR"

done