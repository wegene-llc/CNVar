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
import json


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


def RD_correction(RD, gc_depth_middle):
    RD_correct = []
    RD_correct_rate = []

    for (rd, gc) in RD:
        if gc in gc_depth_middle and gc_depth_middle[gc] != 0:
            rd_correct = rd * gc_depth_middle[gc]
            rd_correct_rate =  gc_depth_middle[gc]
        else:
            rd_correct = rd
            rd_correct_rate = 1
            #print('Warning: cannot find ref depth for GC ' +  str(gc) )
        RD_correct.append(rd_correct)
        RD_correct_rate.append(rd_correct_rate)

    return RD_correct, RD_correct_rate


def RD_with_GC_correct(CHR, start, end, binsize, depth_file, gc_depth_middle, uniq_region, slide):
    RD_correct = []
    # window size = binsize, slide
    gcs = []
    rds = []
    rdc = []
    RD = getRD(CHR,start,end,slide,depth_file, uniq_region, binsize)
    RD_origin = []

    for rd in RD:
        RD_origin.append(rd[0])

    RD_correct, RD_correct_rate = RD_correction(RD,gc_depth_middle)

    return RD_origin, RD_correct, RD_correct_rate


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


def readDepthMethodMatrixMemOptimizeGCmodifyGCcorrect(depth_file, binsize, outfile, slide, region, gcconfig_file):

    binsize = int(binsize)
    slide = int(slide)
    CHR,start,end = parse_region(region)

    with open(gcconfig_file) as f:
        line = f.readline()
        Dict = json.loads(line)
        gc_depth_middle_ratio = {}
        for i in Dict:
            gc_depth_middle_ratio[float(i)] = Dict[i]

    major_cnv = 'N'
    bed_file = 'N'
    Variant = 'N'

    print('''depth_file = %s
        binsize = %s
        outfile = %s
        slide = %s
        start = %s
        end = %s
        CHR = %s
        GCconfig = %s''' % (depth_file,binsize,outfile,slide,start,end,CHR,gcconfig_file) )

    if bed_file != 'N':
        uniq_region = get_uniq_region(bed_file)
    else:
        uniq_region = []

    RD_origin, RD_v, RD_correct_rate = RD_with_GC_correct(CHR, start-50, end, binsize, depth_file,  gc_depth_middle_ratio, uniq_region, slide)

    RD_mean = sum(RD_origin)/len(RD_origin)
    RD_correct_mean = sum(RD_v)/len(RD_v)

    RD_CN = []
    for rd in RD_origin:
        rd_CN = rd/RD_mean*2
        RD_CN.append(rd_CN)

    RD_v_CN = []
    for rd in RD_v:
        rd_v_CN = rd/RD_correct_mean*2
        RD_v_CN.append(rd_v_CN)

    with open(outfile + '.header.txt','w') as h:
        h.write('variant')
        for i in range(0,int((end-start)/slide)-1):
            position = start + slide * i
            if len(uniq_region) > 0:
                if not is_uniq_region(CHR, position, uniq_region, binsize):
                    continue
            h.write('\t' + str(position) )
        h.write('\n')

    with open(outfile+'.origin.txt','w') as w:
        w.write(depth_file)
        for rd in RD_CN:
            w.write('\t%.2f' % rd)
        w.write('\n')

    with open(outfile+'.GCmodified.txt','w') as w:
        w.write(depth_file)
        for rd in RD_v_CN:
            w.write('\t%.2f' % rd)
        w.write('\n')

    with open(outfile+'.GCprofile.txt','w') as w:
        w.write('variant')
        for i in range(0,int((end-start)/slide)-1):
            position = start + slide * i
            if len(uniq_region) > 0:
                if not is_uniq_region(CHR, position, uniq_region, binsize):
                    continue
            w.write('\t' + str(position) )
        w.write('\n')

        w.write(depth_file)
        for rd in RD_correct_rate:
            w.write('\t%.2f' % rd)
        w.write('\n')

    return 0


if __name__ == '__main__':
    if len(sys.argv) < 7:
        print('usage: python %s <depth> <binsize> <outfile> <slide window> <chr:start-end> <GCconfig>'  % (sys.argv[0]))
        exit(1)

    depth_file = sys.argv[1]
    binsize = sys.argv[2]
    outfile = sys.argv[3]
    slide = sys.argv[4]
    region = sys.argv[5]
    gcconfig_file = sys.argv[6]

    readDepthMethodMatrixMemOptimizeGCmodifyGCcorrect(depth_file, binsize, outfile, slide, region, gcconfig_file)