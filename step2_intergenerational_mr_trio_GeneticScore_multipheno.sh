#!/bin/bash

input_path="/ASD/SPARK/WGS/5.mr/output/TTC"
output_path="/ASD/SPARK/WGS/5.mr/output/GeneticScore"
fam_path="/ASD/SPARK/SPARK_phasing_fam/SPARK_WGS_phasing_vcfsample_trio.fam"
PGS_path="/ASD/SSC/WGS/3.PGS/input/GRCh38"

phenotype="trait"
input_directory=$PGS_path/$phenotype
output_directory=$output_path/$phenotype
mkdir -p $output_directory

python intergenerational_mr_trio_gwas_multiPGS.py GeneticScore_multiphenotype \
	-I $input_path/QC_phased_common_rare_all22chr_TTC.vcf.gz \
	-d $input_directory \
	--fam $fam_path \
	-od $output_directory

echo "process ends at:"
date