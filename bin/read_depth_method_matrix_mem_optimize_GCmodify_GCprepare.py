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
import statsmodels.api as sm
import json

def readDepthMethodMatrixMemOptimizeGCmodifyGCprepare(depth_file, binsize, outfile, slide, region):

    binsize = int(binsize)
    slide = int(slide)
    OUTNAME = outfile
    
    CHR,bg_start,bg_end = parse_region(region)

    major_cnv = 'N'
    bed_file = 'N'
    Variant = 'N'

    print('''depth_file = %s
        binsize = %s
        outfile = %s
        slide = %s
        bg_start = %s
        bg_end = %s
        CHR = %s ''' % (depth_file,binsize,outfile,slide,bg_start,bg_end,CHR) )

    if bed_file != 'N':
        uniq_region = get_uniq_region(bed_file)
    else:
        uniq_region = []

    empty = []
    gc_middle, gc_depth_middle, RD_correct_mean = GC_correct_prepare(CHR, bg_start, bg_end, binsize, depth_file, uniq_region, slide, OUTNAME )

    print("gc_middle: %s\nRD_correct_mean: %s" % (gc_middle,RD_correct_mean))

    gc_depth_middle_ratio = {}
    for gc in gc_depth_middle:
        gc_depth_middle_ratio[gc] = gc_depth_middle[gc_middle] / gc_depth_middle[gc]


    gc_depth_middle_ratio_js = json.dumps(gc_depth_middle_ratio)
    with open(outfile,'w') as w:
        w.write(gc_depth_middle_ratio_js)
    
    return 0


def lowess_correction(gc_depth_middle, OUTNAME):
    gc_depth_lowess = {}
    x = []
    y = []
    for i in gc_depth_middle.keys():
            x.append(i)
            y.append(gc_depth_middle[i])

    lowess=sm.nonparametric.lowess

    z=lowess(y,x,frac=0.2)


    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.switch_backend('agg')

    fig = plt.figure(figsize=(40,5),dpi=200)
    plt.plot(x, y,color='r',lw=1)
    plt.plot(z[:,0],z[:,1],color='k',lw=1)
    plt.savefig('%s.lowess_correction.png' % OUTNAME)


    for i in range(len(z[:,0])):
        gc_depth_lowess[z[:,0][i]] =  z[:,1][i]


    return gc_depth_lowess


def safe_open(infile):
    if infile.endswith('.gz'):
        import gzip
        return gzip.open(infile,'rt')
    else:
        return open(infile)


def calucate_RD_in_bin(f,slide):
    count = 0
    rd_bin_sum = 0
    gc_bin_sum = 0
    while count < slide:
        line = f.readline()
        items = line.strip().split()
        chrom = int(items[0])
        pos = int(items[1])
        d = int(items[2])
        b = items[3]
        rd_bin_sum += d
        if b == 'G' or b == 'C':
            gc_bin_sum += 1
        count+=1

    return rd_bin_sum, gc_bin_sum


def getRD(CHR, start,end,slide,depth_file, uniq_region, binsize):
    RD = []
    with safe_open(depth_file) as f:
        for line in f:
            items = line.strip().split()
            chrom = int(items[0])
            pos = int(items[1])

            if chrom == CHR and pos == (start-1):
                # inititate
                # -|-pre|-mid-|-pro-|-
                pre_rd_bin_sum, pre_gc_bin_sum = calucate_RD_in_bin(f,slide)
                mid_rd_bin_sum, mid_gc_bin_sum = calucate_RD_in_bin(f,slide)

                for i in range(start+2*slide,end,slide):
                    rd_sum = 0
                    gc_sum = 0
                    pro_rd_bin_sum, pro_gc_bin_sum = calucate_RD_in_bin(f,slide)
                    if len(uniq_region) > 0:
                        if not is_uniq_region(CHR, i, uniq_region, binsize):
                            continue
                    rd_sum =  pre_rd_bin_sum + mid_rd_bin_sum + pro_rd_bin_sum
                    gc_sum =  pre_gc_bin_sum + mid_gc_bin_sum + pro_gc_bin_sum

                    pre_rd_bin_sum = mid_rd_bin_sum
                    pre_gc_bin_sum = mid_gc_bin_sum

                    mid_rd_bin_sum = pro_rd_bin_sum
                    mid_gc_bin_sum = pro_gc_bin_sum

                    rd = rd_sum / binsize
                    gc = int(gc_sum / binsize * 100)
                    RD.append((rd,gc))

                return RD


def RD_correction(RD, gc_middle, gc_depth_middle):
    RD_correct = []
    RD_correct_rate = []
    for (rd, gc) in RD:
        if gc_middle in gc_depth_middle and gc in gc_depth_middle and gc_depth_middle[gc] != 0:
            rd_correct = rd * gc_depth_middle[gc_middle] / gc_depth_middle[gc]
            rd_correct_rate =  gc_depth_middle[gc_middle] / gc_depth_middle[gc]
        else:
            rd_correct = rd
            rd_correct_rate = 1
            #print('Warning: cannot find ref depth for GC ' +  str(gc) )
        RD_correct.append(rd_correct)
        RD_correct_rate.append(rd_correct_rate)

    return RD_correct, RD_correct_rate


def GC_correct_prepare(CHR, start, end, binsize, depth_file, uniq_region, slide, OUTNAME):
    RD = getRD(CHR, start,end,slide,depth_file, uniq_region, binsize)
    # GC_correct
    RD_gc = []
    gcs = []
    for i in range(len(RD)):
        gcs.append(RD[i][1])
    gcs.sort()
    gc_middle = gcs[int(len(gcs)/2)]

    gc_depth = {}
    for i in range(0,101):
        gc_depth[i] = []
    for i in range(len(RD)):
        gc_depth[RD[i][1]].append(RD[i][0])

    gc_depth_middle = {}
    for i in gc_depth :
        if len(gc_depth[i]) >= 10:           # modify by Enming HE 20190709
            gc_list = sorted(gc_depth[i])
            gc_depth_middle[i] = gc_list[int(len(gc_list)/2 - 1)]

    RD_correct, RD_correct_rate = RD_correction(RD, gc_middle, gc_depth_middle)

    RD_correct_mean = sum(RD_correct)/len(RD_correct)

    gc_depth_middle_lowess = lowess_correction(gc_depth_middle, OUTNAME)

    return gc_middle, gc_depth_middle_lowess, RD_correct_mean


def RD_with_GC_correct(CHR, start, end, binsize, depth_file, gc_middle, gc_depth_middle, uniq_region):
    RD_correct = []
    # window size = binsize, slide
    gcs = []
    rds = []
    rdc = []
    RD = getRD(CHR, start,end,slide,depth_file, uniq_region, binsize)
    RD_correct, RD_correct_rate = RD_correction(RD, gc_middle, gc_depth_middle)

    return RD_correct, RD_correct_rate


def parse_region(region):
    m = re.search('(\d+):(\d+)-(\d+)', region)
    chrom = int(m.group(1))
    start = int(m.group(2))
    end = int(m.group(3))
    return chrom,start,end


def getRegion(TrueSet, Variant):
    Variants = Variant.split('_')
    regions = []
    with open(TrueSet) as f:
        for line in f:
            v = line.strip().split()[3]
            if Variant == 'ALL':
                regions.append((line.strip().split()[1],line.strip().split()[3]))
            else:
                if v in Variants:
                    regions.append((line.strip().split()[1],line.strip().split()[3]))
    return regions


def get_uniq_region(bedfile):
    uniq_region = []
    with open(bedfile) as f:
        for line in f:
            chrom, start, end, length = [ int(x) for x in line.strip().split() ]
            uniq_region.append((chrom, start, end))
    return uniq_region


def is_uniq_region(chrom, pos, uniq_region, binsize):
    check = 0
    for i in uniq_region:
        if chrom == i[0] and pos >= i[1] and pos <= i[2]-binsize:
            check = 1
    return check


if __name__ == '__main__':
    if len(sys.argv) < 6:
        print('usage: python %s <depth> <binsize> <outfile> <slide window> <chr:start-end>'  % (sys.argv[0]))
        exit(1)

    depth_file = sys.argv[1]
    binsize = int(sys.argv[2])
    outfile = sys.argv[3]
    slide = int(sys.argv[4])
    region = sys.argv[5]

    readDepthMethodMatrixMemOptimizeGCmodifyGCprepare(depth_file, binsize, outfile, slide, region)

   