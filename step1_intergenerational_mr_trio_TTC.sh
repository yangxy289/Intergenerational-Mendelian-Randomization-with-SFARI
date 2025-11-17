#!/bin/bash

input_path="/ASD/SPARK/WGS/2.phasing/output/QC_phased_common_rare_vcf"
output_path="/ASD/SPARK/WGS/5.mr/output/TTC"
fam_path="/ASD/SPARK/SPARK_phasing_fam/SPARK_WGS_phasing_vcfsample_trio.fam"

mkdir -p $output_path

for chr in $(seq 1 22); do
       echo "processing TTC: chr$chr"
       python intergenerational_mr_trio_gwas_multiPGS.py TTC \
               -I $input_path/QC_phased_common_rare_chr$chr.vcf.gz \
               --fam $fam_path \
               -w 4000000 > $output_path/QC_phased_common_rare_chr${chr}_TTC.vcf
       echo "done TTC: chr$chr"

       bcftools view -Oz -o $output_path/QC_phased_common_rare_chr${chr}_TTC.vcf.gz $output_path/QC_phased_common_rare_chr${chr}_TTC.vcf
       echo "done bgzip"
       tabix -p vcf $output_path/QC_phased_common_rare_chr${chr}_TTC.vcf.gz
       echo "done tabix"

       rm $output_path/QC_phased_common_rare_chr${chr}_TTC.vcf
done

>$output_path/QC_phased_common_rare_all22chr_TTC.list
for chr in $(seq 1 22); do
        echo "$output_path/QC_phased_common_rare_chr${chr}_TTC.vcf.gz"
        echo "$output_path/QC_phased_common_rare_chr${chr}_TTC.vcf.gz" >> $output_path/QC_phased_common_rare_all22chr_TTC.list
done
echo "done $output_path/QC_phased_common_rare_all22chr_TTC.list"

bcftools concat --threads 50 -Oz -o $output_path/QC_phased_common_rare_all22chr_TTC.vcf.gz -f $output_path/QC_phased_common_rare_all22chr_TTC.list
echo "done bcftools concat"

bcftools index $output_path/QC_phased_common_rare_all22chr_TTC.vcf.gz
echo "done bcftools index"

