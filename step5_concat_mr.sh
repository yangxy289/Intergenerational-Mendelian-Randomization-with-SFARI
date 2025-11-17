#!/bin/bash

input_path="/ASD/SPARK/WGS/5.mr/output/MR"
output_path="/ASD/SPARK/WGS/5.mr/output/MR"

phenotype="trait"
mr_path=$input_path/$phenotype/trio5712_intergenerationalmr
gt_path=$input_path/$phenotype/trio5712_gt_regression
concat_gt_path=$output_path/$phenotype/${phenotype}_trio5712_regression_gt_concat.csv
concat_h14_path=$output_path/$phenotype/${phenotype}_trio5712_intergenerationalMR_h14_concat.csv
concat_all_path=$output_path/$phenotype/${phenotype}_trio_all_concat.csv

n=1
for gt in $gt_path/*_gt_gender_PC.txt; do
    prefix=$(basename "$gt" _trio5712_gt_gender_PC.txt)
    echo "processing concat gt: $gt"

    if [ $n -eq 1 ]; then
        awk -v p="$prefix" -v OFS="," '
            /^ *Feature/ && !header_printed {
                print "Feature","Coef","Stderr","Pvalue","Zscore","CI_lower","CI_upper","OR","OR_CI_lower","OR_CI_upper","phenotype"
                header_printed=1
                next
            }
            /child_genotype_score|father_genotype_score|maternal_genotype_score/ {
                print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,p
            }
        ' "$gt" > "$concat_gt_path"

        n=0
    else
        awk -v p="$prefix" -v OFS="," '
            /child_genotype_score|father_genotype_score|maternal_genotype_score/ {
                print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,p
            }
        ' "$gt" >> "$concat_gt_path"
    fi
done
echo "done concat gt"

n=1
for mr in $mr_path/*_h14_gender_PC.txt; do
        prefix=$(basename $mr _MR_trio_5712_h14_gender_PC.txt)
        echo "processing concat gt: $mr"

        if [ $n -eq 1 ]; then
        awk -v prefix="$prefix" -v OFS="," '
            NR==1 {print "Feature","Coef","Stderr","Pvalue","Zscore","CI_lower","CI_upper","OR","OR_CI_lower","OR_CI_upper","phenotype"}
            NR>1 && ($1=="h1" || $1=="h2" || $1=="h3" || $1=="h4" || $1=="MY" || $1=="FY") {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,prefix}
        ' "$mr" > "$concat_h14_path"
        n=0
    else
        awk -v prefix="$prefix" -v OFS="," '
            NR>1 && ($1=="h1" || $1=="h2" || $1=="h3" || $1=="h4" || $1=="MY" || $1=="FY") {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,prefix}
        ' "$mr" >> "$concat_h14_path"
    fi

done
echo "done concat h14"

awk 'NR==1{print; next} FNR>1' $concat_gt_path $concat_h14_path > $input_path/$phenotype/tmp_merge.csv

cat > $input_path/$phenotype/feature_order.txt <<EOF
child_genotype_score
father_genotype_score
maternal_genotype_score
h1
h2
h3
h4
MY
FY
EOF

awk -F, 'NR==FNR{rank[$1]=++i; next}
         NR>1{print rank[$1]","$0}' OFS=, $input_path/$phenotype/feature_order.txt $input_path/$phenotype/tmp_merge.csv |
sort -t, -k12,12 -k1,1n | \
cut -d, -f2- > $input_path/$phenotype/tmp_sorted.csv

head -n 1 $concat_gt_path > $concat_all_path
cat $input_path/$phenotype/tmp_sorted.csv >> $concat_all_path
rm $input_path/$phenotype/tmp_merge.csv $input_path/$phenotype/tmp_sorted.csv
echo "done concat all"

sed -i 's/maternal_genotype_score/mother_genotype_score/g' "$concat_all_path"

echo "process ends at:"
date
