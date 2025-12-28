"""
Adapted from: Shujia Huang
Date: 2025-11-7
"""
import argparse
import sys
import gzip

from datetime import datetime
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import norm
from scipy import stats

import os

START_TIME = datetime.now()

def get_beta_value_multi(dir_path):
    traits_beta_values = {}

    for fname in os.listdir(dir_path):
        # OmicsPred_hm_GRCh38 olink somascan nightingale metabolon ukb_eur ukb_multi UKB_PPP
        if not (fname.endswith(".txt") or fname.endswith(".gz")):
        #if not (fname.endswith("_GRCh38.txt") or fname.endswith(".gz")):
            continue  

        trait_name = os.path.splitext(fname)[0]  
        beta_dict = {}

        full_path = os.path.join(dir_path, fname)
        open_func = gzip.open if fname.endswith(".gz") else open

        with open_func(full_path, "rt") as IN:
            for line in IN:
                if line.startswith("rsID") or line.startswith("#") or line.startswith("rsid") or line.startswith("CHROM") or line.startswith("SNP") or line.startswith("MarkerName") or line.startswith("CHR") or line.startswith("Chromosome") or line.startswith("chromosome"):
                    continue

                col = line.strip().split()

                if len(col) < 6:
                    continue
                col[1] = col[1] if col[1].startswith("chr") else "chr" + col[1]
                pos = col[1] + ":" + col[2]
                beta_dict[pos] = [col[3].upper(), col[4].upper(), float(col[5])]

        traits_beta_values[trait_name] = beta_dict
        sys.stderr.write("[INFO] Loaded %d variants for trait: %s \n" % (len(beta_dict), trait_name))
        # print(f"[INFO] Loaded {len(beta_dict)} variants for trait: {trait_name}")

    return traits_beta_values


def load_fam_data(fname):
    fam = {}
    with gzip.open(fname, "rt") if fname.endswith(".gz") else open(fname, "rt") as IN:
        for line in IN:
            if line.startswith("#"):
                continue

            col = line.strip().split()
            sid, fid, mid = col[0], col[1], col[2]
            fam[sid] = [sid, fid, mid]

    return fam


def paternal_allele_origin_by_duo(sample_gt, parent_gt, is_paternal_gt=False): 
    """Determine the paternal allele index in `sample_gt` by parent-offspring duos.

    :param sample_gt: The genotype of sample
    :param parent_gt: The genotype of paternal or maternal
    :param is_paternal_gt: The `parent_gt` is from paternal. default: False
    :return:
    """
    # Default value
    is_error_genotype_match = False 
    paternal_allele_origin = [0, False]  # [paternal_genotype_index, is_clear_origin] 

    # Genotype should be: ["0", "0"], ["0", "1"], ["1", "0"] or ["1", "1"]
    s_gt = sample_gt.split("|")
    p_gt = parent_gt.split("|")

    s_gt_sum = sum(map(int, s_gt))  # should be: 0, 1, or 2
    p_gt_sum = sum(map(int, p_gt))  # should be: 0, 1, or 2

    if s_gt_sum == 1 and p_gt_sum == 1: 
        return is_error_genotype_match, paternal_allele_origin

    if p_gt_sum == 0 and s_gt_sum == 0:  # MOM/DAD: 0|0, KID: 0|0 
        paternal_allele_origin = [0, False]

    elif p_gt_sum == 0 and s_gt_sum == 1:  # MOM/DAD: 0|0, KID: 0|1 or 1|0
        if is_paternal_gt:
            paternal_allele_origin = [0, True] if s_gt[0] == "0" else [1, True]
        else:
            # If first allele is the maternal allele, the second one could only be paternal
            paternal_allele_origin = [0, True] if s_gt[0] == "1" else [1, True]

    elif p_gt_sum == 1 and (s_gt_sum == 0 or s_gt_sum == 2):  # MOM/DAD: 0|1 or 1|0, KID: 0|0, 1|1
        paternal_allele_origin = [0, False] 

    elif p_gt_sum == 2 and s_gt_sum == 1:  # MOM/DAD: 1|1, KID: (0|1 or 1|0)
        if is_paternal_gt:
            paternal_allele_origin = [0, True] if s_gt[0] == "1" else [1, True]
        else:
            paternal_allele_origin = [0, True] if s_gt[0] == "0" else [1, True]

    elif p_gt_sum == 2 and s_gt_sum == 2:  # MOM/DAD: 1|1, KID: 1|1
        paternal_allele_origin = [0, False]

    else:
        is_error_genotype_match = True  # probably hit error genotype or de novo muation

    return is_error_genotype_match, paternal_allele_origin


def paternal_allele_origin_by_trio(sample_gt, father_gt, mother_gt): 
    """Determine the paternal allele index in `sample_gt` by parent-offspring duos.

    :param sample_gt: The genotype of sample
    :param father_gt: The genotype of paternal
    :param mother_gt: The genotype of maternal
    """
    # Default value
    is_error_genotype_match = False
    paternal_allele_origin = [0, False]  # [paternal_genotype_index, is_clear_origin]

    # Genotype should be: ["0", "0"], ["0", "1"], ["1", "0"] or ["1", "1"]
    s_gt = sample_gt.split("|")
    f_gt = father_gt.split("|")
    m_gt = mother_gt.split("|")

    s_gt_sum = sum(map(int, s_gt))  # should be: 0, 1, or 2
    f_gt_sum = sum(map(int, f_gt))  # should be: 0, 1, or 2
    m_gt_sum = sum(map(int, m_gt))  # should be: 0, 1, or 2

    if s_gt_sum == 1 and f_gt_sum == 1 and m_gt_sum == 1:  # 0|1, 0|1, 0|1 
        return is_error_genotype_match, paternal_allele_origin

    if s_gt_sum == 0: 
        if f_gt_sum != 2 or m_gt_sum != 2: 
            paternal_allele_origin = [0, False]
        else:
            # DAD: 1|1 or MOM: 1|1 => impossible 
            is_error_genotype_match = True

    elif s_gt_sum == 1:  # KID: 0|1 or 1|0
        if f_gt_sum == 0 and (m_gt_sum == 1 or m_gt_sum == 2):  
            paternal_allele_origin = [0, True] if s_gt[0] == "0" else [1, True]

        elif f_gt_sum == 1 and m_gt_sum == 0:  # DAD: 0|1 or 1|0, MOM: 0|0
            paternal_allele_origin = [0, True] if s_gt[0] == "1" else [1, True]

        elif f_gt_sum == 1 and m_gt_sum == 1:  # DAD: 0|1 or 1|0, MOM: 0|1 or 1|0 
            pass  # has returned the default value

        elif f_gt_sum == 1 and m_gt_sum == 2:  # DAD: 0|1 or 1|0, MOM: 1|1
            paternal_allele_origin = [0, True] if s_gt[0] == "0" else [1, True]

        elif f_gt_sum == 2 and (m_gt_sum == 0 or m_gt_sum == 1):  # DAD: 1|1, MOM: 0|0 or 0|1 or 1|0
            paternal_allele_origin = [0, True] if s_gt[0] == "1" else [1, True]

        else:  # (0|0, 0|0), (1|1, 1|1) 
            is_error_genotype_match = True

    else:  # KID: 1|1 
        if f_gt_sum != 0 and m_gt_sum != 0:  # DAD: 0|1 or 1|1, MOM: 0|1 or 1|1
            paternal_allele_origin = [0, False]

        else:  # DAD == 0|0 or MOM == 0|0
            is_error_genotype_match = True 

    return is_error_genotype_match, paternal_allele_origin


def offspring_genotype_origin(data, fam_idx, index2sample):
    ind_format = {name: i for i, name in enumerate(data[0][8].split(":"))} 
    if "GT" not in ind_format:
        raise ValueError("[ERROR] VCF ERROR: GT is not in FORMAT.") 

    paternal_allele_origin = {}  # key value is the array index of child in VCF 
    for d in data:
        for c, f, m in fam_idx:  
            child = d[c].split(":") 
            father = d[f].split(":") if f is not None else None 
            mother = d[m].split(":") if m is not None else None 

            if (("." in child[ind_format["GT"]]) or
                    (father and "." in father[ind_format["GT"]]) or
                    (mother and "." in mother[ind_format["GT"]])):
                continue

            if (("/" in child[ind_format["GT"]]) or 
                    (father and "/" in father[ind_format["GT"]]) or
                    (mother and "/" in mother[ind_format["GT"]])):
                raise ValueError("[ERORR] Unphased sample: %s, %s or %s. Detail: %s" % (
                    index2sample[f] if f is not None else "",
                    index2sample[m] if m is not None else "",
                    index2sample[c],
                    "\t".join(d[:7] + [d[f] if f is not None else "-",
                                       d[m] if m is not None else "-",
                                       d[c]])))

            # Genotype should be: "0|0", "0|1", "1|0" or "1|1"
            child_gt = child[ind_format["GT"]] 
            father_gt = father[ind_format["GT"]] if father else None
            mother_gt = mother[ind_format["GT"]] if mother else None

            if c not in paternal_allele_origin: 
                # Key value is the index of child in VCF line. [genotype_index, is_clear_origin]
                is_error_genotype_match, paternal_allele_origin[c] = False, [None, False]

            if paternal_allele_origin[c][1]:
                continue

            if father_gt is None or mother_gt is None: 
                if mother_gt: 
                    # maternal-child pair
                    is_error_genotype_match, paternal_allele_origin[c] = paternal_allele_origin_by_duo(
                        child_gt, mother_gt, is_paternal_gt=False)
                elif father_gt:
                    # paternal-child pair
                    is_error_genotype_match, paternal_allele_origin[c] = paternal_allele_origin_by_duo(
                        child_gt, father_gt, is_paternal_gt=True)
                else: 
                    # Single individual
                    is_error_genotype_match, paternal_allele_origin[c] = False, [0, False] 
            else:
                is_error_genotype_match, paternal_allele_origin[c] = paternal_allele_origin_by_trio(
                    child_gt, father_gt, mother_gt)

            if is_error_genotype_match: 
                paternal_allele_origin[c] = [0, False]
                sys.stderr.write("[WARNING] Genotype match failed but still set original and "
                                 "continue: %s \t(father: %s, %s), (mother: %s, %s) and "
                                 "(child: %s, %s).\n" % (
                                     "\t".join(d[0:5]),
                                     index2sample[f] if f is not None else "-",
                                     d[f] if f is not None else "-",
                                     index2sample[m] if m is not None else "-",
                                     d[m] if m is not None else "-",
                                     index2sample[c],
                                     d[c]))
    return paternal_allele_origin


def output_origin_phased(data, paternal_allele_origin):
    ind_format = {name: i for i, name in enumerate(data[0][8].split(":"))}
    for d in data:
        for k, c in paternal_allele_origin.items():
            ind_info = d[k].split(":")  
            if "." in ind_info[ind_format["GT"]]:  # Missing call, do nothing
                continue

            try:
                # adjust the GT to be "Paternal|Maternal"
                gt = ind_info[ind_format["GT"]].split("|")
                ind_info[ind_format["GT"]] = "|".join([gt[c[0]], gt[1 - c[0]]]) 
                d[k] = ":".join(ind_info) 
            except IndexError as e:
                raise ValueError("[ERROR] IndexError: %s\n\n[Target] %s\n[ALL] %s" %
                                 (e, d[k], "\t".join(d)))

        print("%s" % "\t".join(d))

    return


def determine_variant_parent_origin(in_vcf_fn, fam, window=10000):
    """Transform the phased Child genotype from 'Haplotype-Block-A|Haplotype-Block-B'
    to 'Paternal-Haplotype|Maternal-Haplotype'. In a case:

    KID DAD MOM
    0|1 1|1 0|0

    Transform to be:

    KID DAD MOM
    1|0 1|1 0|0

    :param in_vcf_fn:
    :param child_mother_pairs:
    :param window: The size of phasing block in Beagle.
    :return:
    """
    sample2index, index2sample = {}, {} 
    fam_idx = [] 
    n = 0
    data_buffer = []
    prewindow = {} 
    with gzip.open(in_vcf_fn, "rt") if in_vcf_fn.endswith(".gz") else open(in_vcf_fn, "rt") as IN: 
        # VCF file
        for line in IN: 
            if line.startswith("##"): 
                print(line.strip())
                continue

            col = line.strip().split() 
            if line.startswith("#CHROM"): 
                print(line.strip())
                for i in range(9, len(col)):  
                    sample2index[col[i]] = i #sample1:9
                    index2sample[i] = col[i] #9：sample1

                for i in range(9, len(col)):
                    sample_id = col[i]
                    if sample_id in fam: 
                        sid, fid, mid = fam[sample_id]
                        if (fid != "0" and fid not in sample2index) or (mid != "0" and mid not in sample2index):
                            raise ValueError("[ERROR] %s or %s not in VCF" % (fid, mid))

                        fam_idx.append([sample2index[sid],
                                        sample2index[fid] if fid != "0" else None,
                                        sample2index[mid] if mid != "0" else None])
                continue 

            n += 1
            if n % 100000 == 0:
                elapse_time = datetime.now() - START_TIME
                sys.stderr.write("[INFO] Processing %d records done, "
                                 "%d seconds elapsed\n" % (n, elapse_time.seconds))

            if "," in col[4]:  # ignore multi-allelic
                continue

            chrom, pos = col[0], int(col[1])
            window_num = pos // window  # The Phased window or block
            if chrom not in prewindow:
                prewindow[chrom] = window_num

            if len(data_buffer) and (prewindow[chrom] != window_num or data_buffer[0][0] != chrom):
                paternal_allele_origin_idx = offspring_genotype_origin(data_buffer, fam_idx, index2sample)
                output_origin_phased(data_buffer, paternal_allele_origin_idx)

                data_buffer = []  # clear record
                prewindow[chrom] = window_num  # reset window

            data_buffer.append(col)  # each row is one line of VCF record

        if len(data_buffer):
            paternal_allele_origin_idx = offspring_genotype_origin(data_buffer, fam_idx, index2sample)
            output_origin_phased(data_buffer, paternal_allele_origin_idx)

    elapse_time = datetime.now() - START_TIME
    sys.stderr.write("[INFO] All %d records loaded, %d seconds elapsed.\n" % (n, elapse_time.seconds))

    return


def calculate_genotype_and_haplotype_score_multi(in_vcf_fn, traits_beta_values, fam, outdir, score_model, is_dosage=False):
    sample2index, index2sample = {}, {}
    gs = {}  # { sample_key: { trait: [fet_sum, fat_sum, mat_sum, h1_sum, h2_sum, h3_sum, h4_sum, count] } }
    child_father_mother_idx = []
    n = 0
    with gzip.open(in_vcf_fn, "rt") if in_vcf_fn.endswith(".gz") else open(in_vcf_fn, "rt") as IN:
        # VCF file
        for line in IN:
            if line.startswith("##"):
                continue

            col = line.strip().split()
            if line.startswith("#CHROM"):
                for i in range(9, len(col)):  # load sample ID and the index of sample
                    sample2index[col[i]] = i
                    index2sample[i] = col[i]

                for i in range(9, len(col)): 
                    sample_id = col[i]
                    if sample_id in fam:
                        sid, fid, mid = fam[sample_id]
                        if mid == "0" and fid == "0":
                           continue
                        elif fid != "0" and mid == "0":
                           child_father_mother_idx.append([sample2index[sid], sample2index[fid], "0"])
                        elif fid == "0" and mid != "0":
                           child_father_mother_idx.append([sample2index[sid], "0", sample2index[mid]])
                        else:
                           child_father_mother_idx.append([sample2index[sid], sample2index[fid], sample2index[mid]])

                continue

            n += 1
            if n % 100000 == 0:
                elapse_time = datetime.now() - START_TIME
                sys.stderr.write("[INFO] Processing %d records done, %d seconds elapsed\n" % (n, elapse_time.seconds))

            pos = col[0] + ":" + col[1]
            ref_allele = col[3].upper()
            alt_allele = col[4].upper()

            if "," in alt_allele:
                continue

            ind_format = {name: i for i, name in enumerate(col[8].split(":"))}
            if "GT" not in ind_format:
                raise ValueError("[ERROR] 'GT' filed is required in VCF for each individual.")

            if is_dosage and ("GP" not in ind_format) and ("DS" not in ind_format):
                raise ValueError("[ERROR] 'GP' or 'DS' field is required for dosage "
                                 "for each individual.")

            for trait, beta_dict in traits_beta_values.items():
                if pos not in beta_dict:
                    continue

                a1, a2, beta = beta_dict[pos]
                if (a2 != "-") and (ref_allele + alt_allele != a1 + a2) and \
                        (alt_allele + ref_allele != a1 + a2):
                    #sys.stderr.write("[ERROR] Alleles not matched: [%s, %s] != [%s, %s]\n"
                    #                 "%s" % (ref_allele, alt_allele, a1, a2, line))
                    sys.stderr.write("[ERROR] Alleles not matched: [%s, %s] != [%s, %s]\n"
                                     "%s\t" "%s\t" "%s\t" "%s\t" "%s\n" % (ref_allele, alt_allele, a1, a2, col[0], col[1], col[2], col[3], col[4]))

                    continue
                if ref_allele == a1:
                    beta = -beta

                for c, f, m in child_father_mother_idx:
                    if f != "0" and m == "0":
                        father_gt_str = col[f].split(":")[ind_format["GT"]]
                        child_gt_str = col[c].split(":")[ind_format["GT"]]
                        if ("." in father_gt_str) or ("." in child_gt_str): continue
                        father_gt = list(map(int, father_gt_str.split("|")))
                        child_gt = list(map(int, child_gt_str.split("|")))
                        if (sum(father_gt) == 0 and sum(child_gt) == 2) or (
                                sum(father_gt) == 2 and sum(child_gt) == 0): continue
                        fet = sum(child_gt)
                        fat = sum(father_gt)
                        mat = 0
                        h1 = child_gt[1]
                        h2 = 0
                        h3 = child_gt[0]
                        h4 = father_gt[0] if father_gt[0] != child_gt[0] else father_gt[1]
                        k = index2sample[c] + "_" + index2sample[f] + "_0"

                    elif f == "0" and m != "0":
                        mother_gt_str = col[m].split(":")[ind_format["GT"]]
                        child_gt_str = col[c].split(":")[ind_format["GT"]]
                        if ("." in mother_gt_str) or ("." in child_gt_str): continue
                        mother_gt = list(map(int, mother_gt_str.split("|")))
                        child_gt = list(map(int, child_gt_str.split("|")))
                        if (sum(mother_gt) == 0 and sum(child_gt) == 2) or (
                                sum(mother_gt) == 2 and sum(child_gt) == 0): continue
                        fet = sum(child_gt)
                        fat = 0
                        mat = sum(mother_gt)
                        h1 = child_gt[1]
                        h2 = mother_gt[0] if mother_gt[0] != child_gt[1] else mother_gt[1]
                        h3 = child_gt[0]
                        h4 = 0
                        k = index2sample[c] + "_0_" + index2sample[m]

                    else:
                        father_gt_str = col[f].split(":")[ind_format["GT"]]
                        mother_gt_str = col[m].split(":")[ind_format["GT"]]
                        child_gt_str = col[c].split(":")[ind_format["GT"]]
                        if ("." in father_gt_str) or ("." in mother_gt_str) or ("." in child_gt_str): continue
                        father_gt = list(map(int, father_gt_str.split("|")))
                        mother_gt = list(map(int, mother_gt_str.split("|")))
                        child_gt = list(map(int, child_gt_str.split("|")))
                        # if (sum(father_gt) == 0 and sum(child_gt) == 2) or (sum(father_gt) == 2 and sum(child_gt) == 0) or (sum(mother_gt) == 0 and sum(child_gt) == 2) or (sum(mother_gt) == 2 and sum(child_gt) == 0): continue
                        # if (sum(father_gt) == 0 and sum(child_gt) == 2) or (sum(father_gt) == 2 and sum(child_gt) == 0) or (sum(mother_gt) == 0 and sum(child_gt) == 2) or (sum(mother_gt) == 2 and sum(child_gt) == 0) or (sum(father_gt) == 0 and sum(mother_gt) == 0 and sum(child_gt) == 1) or (sum(father_gt) == 2 and sum(mother_gt) == 2 and sum(child_gt) == 1): continue
                        if (sum(father_gt) == 0 and sum(child_gt) == 2) or (
                                sum(father_gt) == 2 and sum(child_gt) == 0) or (
                                sum(mother_gt) == 0 and sum(child_gt) == 2) or (
                                sum(mother_gt) == 2 and sum(child_gt) == 0) or (
                                sum(father_gt) == 0 and sum(mother_gt) == 0 and sum(child_gt) == 1) or (
                                sum(father_gt) == 2 and sum(mother_gt) == 2 and sum(child_gt) == 1) or (
                                sum(father_gt) == 1 and sum(mother_gt) == 1 and sum(child_gt) == 1): continue
                        fet = sum(child_gt)
                        fat = sum(father_gt)
                        mat = sum(mother_gt)
                        h1 = child_gt[1]
                        h2 = mother_gt[0] if mother_gt[0] != child_gt[1] else mother_gt[1]
                        h3 = child_gt[0]
                        h4 = father_gt[0] if father_gt[0] != child_gt[0] else father_gt[1]
                        k = index2sample[c] + "_" + index2sample[f] + "_" + index2sample[m]

                    s_fet = fet * beta
                    s_fat = fat * beta
                    s_mat = mat * beta
                    s_h1 = h1 * beta
                    s_h2 = h2 * beta
                    s_h3 = h3 * beta
                    s_h4 = h4 * beta
                    if k not in gs: gs[k] = {}
                    if trait not in gs[k]: gs[k][trait] = [0, 0, 0, 0, 0, 0, 0, 0]
                    gs[k][trait][0] += s_fet
                    gs[k][trait][1] += s_fat
                    gs[k][trait][2] += s_mat
                    gs[k][trait][3] += s_h1
                    gs[k][trait][4] += s_h2
                    gs[k][trait][5] += s_h3
                    gs[k][trait][6] += s_h4
                    gs[k][trait][7] += 1

    elapse_time = datetime.now() - START_TIME
    sys.stderr.write("[INFO] All %d records loaded, %d seconds elapsed.\n" % (n, elapse_time.seconds))

    os.makedirs(outdir, exist_ok=True)

    trait_results = {}
    for k, traits in gs.items():
        for trait, values in traits.items():
            if trait not in trait_results: trait_results[trait] = {}
            trait_results[trait][k] = values

    for trait, k in trait_results.items():
        outfile = os.path.join(outdir, f"{trait}.txt")
        with open(outfile, "w") as f:
            f.write("#Child\tFather\tMother\tchild_genotype_score\tfather_genotype_score\tmaternal_genotype_score\th1\th2\th3\th4\tsite_number\n")
            for sample_key, values in k.items():
                child, father, mother = sample_key.split("_", 2)
                if score_model == "avg" and values[7] > 0: score = [x / values[7] for x in values[:-1]]
                else: score = values[:-1]
                score_str = "\t".join(map(str, score))
                f.write(f"{child}\t{father}\t{mother}\t{score_str}\t{values[7]}\n")

        sys.stderr.write("[INFO] Trait %s results saved to %s \n" % (trait, outfile))

    return


def phenotype_concat(in_gs_fn, in_pheno_file): 
    gs_data = {}
    header = []
    with open(in_gs_fn, "rt") as I:
        """
        #Sample_G1	Sample_G2	maternal_genotype_score	child_genotype_score	M1(C1)	M2	C2
        00113051204M47BFF2	13120631BFF2	-18.16394	-18.3755	-10.6309	0.000978	-10.5710
        00113061190M24BFF2	13121060BFF2	-18.17056	-18.2710	-10.6823	-0.00194	-10.5112
        """
        n = 0
        for line in I:
            n += 1
            if n == 1:
                header = line.strip().split()
                continue

            col = line.strip().split()
            gs_data[col[0]] = col

    with open(in_pheno_file, "rt") as I:
        n = 0
        for line in I:
            """
            IID	AGE	SEX	height	weight	PREBMI	Weight_gain	......
            00113051204M47BFF2	29	2	159.0	44.0	17	......
            """
            n += 1
            if n == 1:
                header += line.strip().split()[1:]
                print("%s" % "\t".join(header))
                continue

            col = line.strip().split()

            sample_id = col[0]
            if sample_id in gs_data:
                print("%s" % "\t".join(gs_data[sample_id] + col[1:]))


def intergenerationalMR_regression(data, y_name, x_names, covar_names=""):
    d = x_names.strip().split(",")
    c = covar_names.strip().split(",") if covar_names.strip() else []
    X = data[d + c].copy()  
    y = data[y_name]

    for col in X.columns:
        if pd.api.types.is_numeric_dtype(X[col]) and X[col].nunique() > 2:
            X[col] = (X[col] - X[col].mean()) / X[col].std()
        else:
            X[col] = X[col]
    X = sm.add_constant(X)  

    if pd.api.types.is_numeric_dtype(y):
        unique_y = np.sort(y.unique())
        if len(unique_y) == 2:
            label_mapping = {unique_y[0]: 0, unique_y[1]: 1}
            y = y.map(label_mapping)
            regression = sm.Logit(y, X)
            model_type = f"Logistic (labels mapped: {unique_y[0]}->0, {unique_y[1]}->1)"
        elif len(unique_y) > 2 and all(y.dropna().astype(int) == y.dropna()):
            regression = sm.MNLogit(y, X)
            model_type = "Multinomial Logistic"
        else:
            regression = sm.OLS(y, X)
            model_type = "OLS"
    else:
        raise ValueError("Not numeric")

    model = regression.fit()

    conf_int = model.conf_int()
    conf_int.columns = ["CI_lower", "CI_upper"]

    fe = pd.concat(
        [model.params, model.bse, model.pvalues, model.tvalues, conf_int],
        axis=1
    ).reset_index().rename(
        columns={
            "index": "Feature",
            0: "Coef",
            1: "Stderr",
            2: "Pvalue",
            3: "Zscore"
        }
    ).sort_values(by=["Pvalue"], inplace=False)
    fe = fe[fe["Feature"].isin(["h1", "h2", "h3", "h4"])].reset_index(drop=True)

    # MY、FY
    m_beta = 0.5 * (float(fe[fe["Feature"] == "h1"].iloc[0]["Coef"]) +
                    float(fe[fe["Feature"] == "h2"].iloc[0]["Coef"]) -
                    float(fe[fe["Feature"] == "h3"].iloc[0]["Coef"]))
    m_se = 0.5 * np.sqrt(float(fe[fe["Feature"] == "h1"].iloc[0]["Stderr"]) ** 2 +
                         float(fe[fe["Feature"] == "h2"].iloc[0]["Stderr"]) ** 2 +
                         float(fe[fe["Feature"] == "h3"].iloc[0]["Stderr"]) ** 2)
    m_t = m_beta / m_se
    m_p = 2 * (1 - stats.t.cdf(abs(m_t), len(X) - 1))   # t
    f_beta = 0.5 * (float(fe[fe["Feature"] == "h1"].iloc[0]["Coef"]) -
                    float(fe[fe["Feature"] == "h2"].iloc[0]["Coef"]) +
                    float(fe[fe["Feature"] == "h3"].iloc[0]["Coef"]))
    f_se = 0.5 * np.sqrt(float(fe[fe["Feature"] == "h1"].iloc[0]["Stderr"]) ** 2 +
                         float(fe[fe["Feature"] == "h2"].iloc[0]["Stderr"]) ** 2 +
                         float(fe[fe["Feature"] == "h3"].iloc[0]["Stderr"]) ** 2)
    f_t = f_beta / f_se
    f_p = 2 * (1 - stats.t.cdf(abs(f_t), len(X) - 1))   # t
    fe = pd.concat([
        fe,
        pd.DataFrame([{"Feature": "MY", "Coef": m_beta, "Stderr": m_se, "Pvalue": m_p, "Zscore": m_t, "CI_lower": m_beta - 1.96*m_se, "CI_upper": m_beta + 1.96*m_se}]),
        pd.DataFrame([{"Feature": "FY", "Coef": f_beta, "Stderr": f_se, "Pvalue": f_p, "Zscore": f_t, "CI_lower": f_beta - 1.96*f_se, "CI_upper": f_beta + 1.96*f_se}])
    ], ignore_index=True)

    # OR
    if model_type.startswith("Logistic"):
        fe["OR"] = np.exp(fe["Coef"])
        fe["OR_CI_lower"] = np.exp(fe["Coef"] - 1.96 * fe["Stderr"])
        fe["OR_CI_upper"] = np.exp(fe["Coef"] + 1.96 * fe["Stderr"])
    else:
        fe["OR"] = np.nan
        fe["OR_CI_lower"] = np.nan
        fe["OR_CI_upper"] = np.nan

    fe = fe.round(6)
    print(fe.to_string(index=False))

    return

def genotype_asd_regression(data, y_name, x_names, covar_names=""):
    d = x_names.strip().split(",")
    c = covar_names.strip().split(",") if covar_names.strip() else []
    X = data[d + c].copy()  
    y = data[y_name]

    for col in X.columns:
        if pd.api.types.is_numeric_dtype(X[col]) and X[col].nunique() > 2:
            X[col] = (X[col] - X[col].mean()) / X[col].std()
        else:
            X[col] = X[col]
    X = sm.add_constant(X)

    if pd.api.types.is_numeric_dtype(y):
        unique_y = np.sort(y.unique())
        if len(unique_y) == 2:
            label_mapping = {unique_y[0]: 0, unique_y[1]: 1}
            y = y.map(label_mapping)
            regression = sm.Logit(y, X)
            model_type = f"Logistic (labels mapped: {unique_y[0]}->0, {unique_y[1]}->1)"
        elif len(unique_y) > 2 and all(y.dropna().astype(int) == y.dropna()):
            regression = sm.MNLogit(y, X)
            model_type = "Multinomial Logistic"
        else:
            regression = sm.OLS(y, X)
            model_type = "OLS"
    else:
        raise ValueError("Not numeric")

    model = regression.fit()

    conf_int = model.conf_int()
    conf_int.columns = ["CI_lower", "CI_upper"]
    fe = pd.concat([model.params, model.bse, model.pvalues, model.tvalues, conf_int], axis=1).reset_index().rename(columns={"index": "Feature", 0: "Coef", 1: "Stderr", 2: "Pvalue", 3: "Zscore"}).sort_values(by=["Pvalue"], inplace=False)
    fe = fe[fe["Feature"].isin(["child_genotype_score", "father_genotype_score", "maternal_genotype_score"])].reset_index(drop=True)

    if model_type.startswith("Logistic"):
        fe["OR"] = np.exp(fe["Coef"])
        fe["OR_CI_lower"] = np.exp(fe["Coef"] - 1.96 * fe["Stderr"])
        fe["OR_CI_upper"] = np.exp(fe["Coef"] + 1.96 * fe["Stderr"])
    else:
        fe["OR"] = np.nan
        fe["OR_CI_lower"] = np.nan
        fe["OR_CI_upper"] = np.nan

    fe = fe.round(6)
    print(fe.to_string(index=False))

    return


if __name__ == "__main__":
    cmd_parser = argparse.ArgumentParser(description="Usage: ")
    commands = cmd_parser.add_subparsers(dest="command", title="Commands")
    tc_cmd = commands.add_parser("TTC", help="Transform the phased genotype of child from "
                                             "'Haplotype-Block-A|Haplotype-Block-B' to "
                                             "'Paternal-Haplotype|Maternal-Haplotype'.")
    tc_cmd.add_argument("-I", "--target", dest="target", type=str, required=True,
                        help="Input a phased VCF. Required.")
    tc_cmd.add_argument("-w", "--window", dest="window", type=int, default=10000, required=False,
                        help="The phased block size. [10000]")
    tc_cmd.add_argument("--fam", dest="fam", type=str, required=True,
                        help="Input a .fam file with mother and children.")

    gs_cmd = commands.add_parser("GeneticScore_multiphenotype", help="Calculate multiple genotypes Genetic score.")
    gs_cmd.add_argument("-I", "--target", dest="target", type=str, required=True,
                        help="Input VCF. Required.")
    gs_cmd.add_argument("-d", "--indirectory", dest="indirectory", type=str, required=True,
                        help="A directory with GWAS pos beta files of phenotypes")
    gs_cmd.add_argument("--fam", dest="fam", type=str, required=False,
                        help="Input a .fam file [option]. If provide .fam file, this module will only "
                             "calculate the genetic score for mother-child pairs according to the "
                             "parent_origin VCF, which create by 'TTC' module.")
    gs_cmd.add_argument("-od", "--outdirectory", dest="outdirectory", type=str, required=True,
                        help="A directory with geneticscore files of phenotypes")
    gs_cmd.add_argument("--score-model", dest="sm", type=str, default="avg", required=False,
                        help="The model for calculating genetic score [avg, sum]. Default: avg.")
    gs_cmd.add_argument("--dosage", dest="dosage", action="store_true", help="Use dosage.")

    mr_cmd = commands.add_parser("MR", help="Mendelian Randomization")
    mr_cmd.add_argument("-I", "--input", dest="input", type=str, required=True,
                        help="Input data file")
    mr_cmd.add_argument("-x", dest="x_name", type=str, required=True,
                        help="Load the designated phenotype(s) as x from the '--input'.")
    mr_cmd.add_argument("--covar", dest="covar_name", type=str, required=False, default="", 
                        help="Only load the designated covariate(s) from the '--input'. "
                        "If not set, no covariates will be used.")
    mr_cmd.add_argument("-y", dest="y_name", type=str, required=True,
                        help="Load the designated data as y from the '--input'.")

    add_cmd = commands.add_parser("ADD", help="Concat genetic score data together with phenotype data.")
    add_cmd.add_argument("-g", dest="genetic_score", type=str, required=True, help="input genetic score file.")
    add_cmd.add_argument("-p", dest="phenotype", type=str, required=True, help="input phenotype data file.")

    args = cmd_parser.parse_args()
    if args.command == "TTC":
        fam_data = load_fam_data(args.fam) 
        determine_variant_parent_origin(args.target, fam_data, window=args.window)

    elif args.command == "GeneticScore_multiphenotype":

        if args.sm not in ["avg", "sum"]:
            raise ValueError("[ERROR] The value of argument '--score-model' could only be 'avg' or 'sum'.")

        if args.fam:
            fam_data = load_fam_data(args.fam)
            traits_beta_values = get_beta_value_multi(args.indirectory)
            calculate_genotype_and_haplotype_score_multi(args.target, traits_beta_values, fam_data, args.outdirectory, score_model=args.sm, is_dosage=args.dosage)
        else:
            beta_value = get_beta_value(args.base)
            calculate_genotype_score(args.target, beta_value, score_model=args.sm, is_dosage=args.dosage)

    elif args.command == "MR":
        data = pd.read_table(args.input, sep="\t")
        intergenerationalMR_regression(data, args.y_name, args.x_name, args.covar_name)
        #genotype_asd_regression(data, args.y_name, args.x_name, args.covar_name)

    elif args.command == "ADD":
        phenotype_concat(args.genetic_score, args.phenotype)

    else:
        cmd_parser.print_help()

    elapsed_time = datetime.now() - START_TIME
    sys.stderr.write("\n** process done, %d seconds elapsed **\n" % elapsed_time.seconds)



