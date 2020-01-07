#! /usr/bin/env python
# -*- coding: utf-8 -*-
#########################################################################################################
#    CNVar: accurate knowledge-based copy number variation prediction using max likelihood method       #
#                                                                                                       # 
#    Copyright (C) 2019 Enming He(emhe@wegene.com)                                                      #
#                                                                                                       #
#    This program is free software: you can redistribute it and/or modify                               #
#    it under the terms of the GNU General Public License as published by                               #
#    the Free Software Foundation, either version 3 of the License, or                                  #
#    (at your option) any later version.                                                                #
#                                                                                                       #
#    This program is distributed in the hope that it will be useful,                                    #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of                                     #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                      #
#    GNU General Public License for more details.                                                       #
#                                                                                                       #
#    You should have received a copy of the GNU General Public License                                  #
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.                             #
#########################################################################################################

import os
import sys
import re
import numpy as np
import scipy.stats as stats
import math
import subprocess
import argparse
import read_depth_method_matrix_mem_optimize_GCmodify_GCcorrect


def safe_open(infile):
    if infile.endswith('.gz'):
        import gzip
        return gzip.open(infile)
    else:
        return open(infile)


def read_GCprofile(infile):
    print('reading:GCprofile %s ...' % infile)
    GCprofile = {}
    with open(infile) as f:
        line = f.readline()
        headers = line.strip().split()
        line = f.readline()
        gc_values = line.strip().split()
    
    for i in range(1,len(headers)):
        GCprofile[int(headers[i])] = float(gc_values[i])
    print('finish reading:GCprofile')
    return GCprofile


def possion(cov):
    if cov <= 60:
        rate = 60
    else:
        rate = int(cov)
    Dict = {}
    n = np.arange((rate+1)*2)
    y = stats.poisson.pmf(n,cov)

    for i in n:
        Dict[i] = y[i]
    return Dict


def reverse_possion(cov):
    rate = int(cov)
    Dict = {}
    n = np.arange(rate)
    y = stats.poisson.pmf(n,rate)
    y1 = y.tolist()
    y1.reverse()
    for i in n:
        Dict[i] = y1[i]
    return Dict


def log(Dict,x,logmin):
    if x > max(Dict.keys()) or x < min(Dict.keys()):
        return math.log(logmin)
    else:
        if Dict[x] < logmin:
            return math.log(logmin)
        else:
            return math.log(Dict[x])


def possibility(D,DR,cov,N, CHROM_TEST,START_TEST,END_TEST, MISMATCH, logmin):
    lgP = math.log(1)
    CN_models = {}
    
    def possionCN(CN):
        if (CN!=0):
            CN_model = possion(cov*CN/2)
        else:
            CN_model  = possion(MISMATCH * cov)
        return CN_model
    
    for i in range(START_TEST,END_TEST):
        depth = D[(CHROM_TEST,i)]
        CN = round(DR[(CHROM_TEST,i)]/N,2)
        
        if CN not in CN_models:
            CN_models[CN] = possionCN(CN)
    
        y = log(CN_models[CN],depth, logmin)
       
        lgP += y

    return lgP


def read_ref(reffile, CHROM_TEST, START_TEST,END_TEST):
    print('reading ref file: %s ' % reffile)
    depth_sum =0
    count = 0
    D = {}
    with safe_open(reffile) as f:
        for line in f:
            chrom, pos, depth = [int(x) for x in line.strip().split()[0:3]]
            if chrom == CHROM_TEST and pos<= END_TEST and pos>= START_TEST:
                D[(chrom,pos)] = depth
            elif chrom == CHROM_TEST and pos > END_TEST:
                break
    return D


def read_sample(infile,GCprofile, CHROM_TEST, START_TEST,END_TEST ):
    print('reading sample file: %s ...' % infile)
    with safe_open(infile) as f:
        D = {}
        depth_sum = 0
        count = 0
        for line in f:
            chrom, pos, depth_raw = [int(x) for x in line.strip().split()[0:3]]
            if chrom == CHROM_TEST and pos<= END_TEST and pos>= START_TEST:
                depth = int(depth_raw * GCprofile[pos - pos % 50])
                D[(chrom,pos)] = depth
            elif chrom == CHROM_TEST and pos > END_TEST:   # if assume background region behinds testing region
                break
    return D


def calculator(LIST,sample,reflist, binsize, slide, gc_config,CHROM_TEST, START_TEST, END_TEST, AVERAGE_DEPTH, MISMATCH, logmin):
    min_start =  START_TEST
    max_end = END_TEST
    
    read_depth_method_matrix_mem_optimize_GCmodify_GCcorrect.readDepthMethodMatrixMemOptimizeGCmodifyGCcorrect(sample, binsize, sample, slide, '%s:%s-%s' % (CHROM_TEST, START_TEST-150, END_TEST+150), gc_config) # 150 , end should be wilder

    GCprofile = read_GCprofile('%s.GCprofile.txt' % sample)
    D = read_sample(sample,GCprofile, CHROM_TEST, START_TEST, END_TEST)
    cov = AVERAGE_DEPTH
    for ref_info in reflist:
        ref, Depth_Ref = ref_info.strip().split()
        N = float(Depth_Ref)/2
        DR = read_ref(ref,CHROM_TEST, START_TEST,END_TEST)
        lgP = possibility(D, DR, cov, N, CHROM_TEST,START_TEST,END_TEST, MISMATCH,logmin)
        LIST.append((ref,lgP))


def possionFixPossionModelSiteBYSiteOptionGCmodify(reference_list, sample_list, outfile, threads, detect_region, average_depth, gc_config, binsize, slide_window, logmin, mismatch ):
    
    logmin = float(logmin)
    m1 = re.search('(\d+):(\d+)-(\d+)',detect_region)
    CHROM_TEST = int(m1.group(1))
    START_TEST = int(m1.group(2))
    END_TEST = int(m1.group(3))
    AVERAGE_DEPTH = float(average_depth)
    MISMATCH = float(mismatch)
    
    print('''
        CHROM_TEST = %s
        START_TEST = %s
        END_TEST = %s
        ''' % (CHROM_TEST, START_TEST , END_TEST)
          )
        
    reflist=[]
    samplelist=[]
    
    with open(reference_list) as f:
        for line in f:
            reflist.append(line.strip())

    with open(sample_list) as f:
        for line in f:
            samplelist.append(line.strip())
        
        sample_Dict = {}
        
        for sample in samplelist:
            sample_Dict[sample] = []
            calculator(sample_Dict[sample],sample,reflist, binsize, slide_window,
            gc_config,CHROM_TEST, START_TEST, END_TEST, AVERAGE_DEPTH, MISMATCH, logmin)

    with open(outfile,'w') as w:
        for sample in sample_Dict:
            mylist = sorted(sample_Dict[sample],key=lambda x: x[1],reverse=True)
            rank = 1
            for i in mylist:
                w.write('%s\t%s\t%s\t%s\n' % (sample,i[0], i[1],rank))
                rank += 1

    return 0


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-rl', '--reference_list', help='reference list',required=True)
    parser.add_argument('-sl', '--sample_list', help='sample list',required=True)
    parser.add_argument('-o', '--outfile', help='outfile',required=True)
    parser.add_argument('-t', '--threads', help='threads',default=16)
    parser.add_argument('-r', '--detect_region', help='detect region, chr:start-end ',required=True)
    parser.add_argument('-d', '--average_depth', help='average depth',required=True)
    parser.add_argument('-c', '--gc_config', help='GC config',required=True)
    parser.add_argument('-bz', '--binsize', help='binsize',required=True)
    parser.add_argument('-sw', '--slide_window', help='slide_window',required=True)
    parser.add_argument('-l', '--logmin', help='logmin', default=1e-7)
    parser.add_argument('-m', '--mismatch', help='assume mismatch', default=1e-4)
    args = parser.parse_args()
    
    reference_list = args.reference_list
    sample_list = args.sample_list
    outfile = args.outfile
    threads = args.threads
    detect_region = args.args.detect_region
    average_depth = args.average_depth
    gc_config = args.gc_config
    binsize = args.binsize
    slide_window = args.slide_window
    logmin = args.logmin
    mismatch =  args.mismatch
    
    possionFixPossionModelSiteBYSiteOptionGCmodify(reference_list, sample_list, outfile, threads, detect_region, average_depth, gc_config, binsize, slide_window, logmin, mismatch )

