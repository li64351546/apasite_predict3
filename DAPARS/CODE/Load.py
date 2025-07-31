import numpy as np
import os
import sys
import datetime

import scipy as sp
import scipy.stats

from bisect import bisect

import math



def Load_Target_Wig_files(All_Sample_files, Annotated_3UTR_file):

    UTR_events_dict = {}
    All_Samples_Total_depth = []

    # Annotated_3UTR_file->UTR_events_dict
    # chr14	50792327	50792946	NM_001003805|ATP5S|chr14|+	0->NM_001003805|ATP5S|chr14|+
    # Annotated_3UTR_file=hg19_refseq_extracted_3UTR.bed 25294条
    for line in open(Annotated_3UTR_file,'r'):
        fields = line.strip('\n').split('\t')
        curr_chr = fields[0]
        region_start = int(float(fields[1]))
        region_end   = int(float(fields[2]))
        curr_strand  = fields[-1]

        UTR_pos = "%s:%s-%s" % (curr_chr, region_start, region_end)
        end_shift = int(round(abs(int(region_start) - int(region_end)) * 0.2))
        if curr_strand == '+':
            region_end = str(int(region_end) - end_shift)
        else:
            region_start = str(int(region_start) + end_shift)
        region_start = int(region_start) + 1
        region_end   = int(region_end) - 1
        if region_start + 50 < region_end:
            UTR_events_dict[fields[3]] = [fields[0],region_start,region_end,fields[-1],UTR_pos]





    ## 加载所有样品的read coverage
    ## Load coverage for all samples
    All_samples_extracted_3UTR_coverage_dict = {}
    for curr_wig_file in All_Sample_files:

        curr_sample_All_chroms_coverage_dict = {}
        num_line = 0
        cur_sample_total_depth = 0



        for line in open(curr_wig_file,'r'):

            if '#' not in line and line[0:3] == 'chr':
                fields = line.strip('\n').split('\t')
                chrom_name = fields[0]
                region_start = int(float(fields[1]))
                region_end = int(float(fields[2]))
                cur_sample_total_depth += int(float(fields[-1])) * (region_end - region_start)
                if chrom_name not in curr_sample_All_chroms_coverage_dict:
                    curr_sample_All_chroms_coverage_dict[chrom_name] = [[0],[0]]

                if region_start > curr_sample_All_chroms_coverage_dict[chrom_name][0][-1]:
                    curr_sample_All_chroms_coverage_dict[chrom_name][0].append(region_start)
                    curr_sample_All_chroms_coverage_dict[chrom_name][1].append(0)
                curr_sample_All_chroms_coverage_dict[chrom_name][0].append(region_end)
                curr_sample_All_chroms_coverage_dict[chrom_name][1].append(int(float(fields[-1])))
            num_line += 1
        curr_sample_All_chroms_coverage_dict[chrom_name][1].append(0)
        All_Samples_Total_depth.append(cur_sample_total_depth)






        for curr_3UTR_event_id in UTR_events_dict:
            curr_3UTR_structure = UTR_events_dict[curr_3UTR_event_id]
            curr_chr = curr_3UTR_structure[0]

            if curr_chr in curr_sample_All_chroms_coverage_dict:
                curr_chr_coverage = curr_sample_All_chroms_coverage_dict[curr_chr]
                region_start = curr_3UTR_structure[1]
                region_end = curr_3UTR_structure[2]
                left_region_index = bisect(curr_chr_coverage[0],region_start)
                right_region_index = bisect(curr_chr_coverage[0],region_end)

                extracted_coverage = curr_chr_coverage[1][left_region_index:right_region_index+1]
                extracted_3UTR_region = curr_chr_coverage[0][left_region_index:right_region_index]
                extracted_3UTR_region.insert(0,region_start)
                extracted_3UTR_region.append(region_end)
                
                if curr_3UTR_event_id not in All_samples_extracted_3UTR_coverage_dict:
                    All_samples_extracted_3UTR_coverage_dict[curr_3UTR_event_id] = []
                All_samples_extracted_3UTR_coverage_dict[curr_3UTR_event_id].append([extracted_coverage,extracted_3UTR_region])



    return All_samples_extracted_3UTR_coverage_dict,np.array(All_Samples_Total_depth),UTR_events_dict