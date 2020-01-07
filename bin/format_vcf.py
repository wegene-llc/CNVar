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

import sys
import re

def cnvar2vcf( outputfile, vcfname, mutation_info ):
    Dict = {}
    with open(mutation_info) as f:
        for line in f:
            if line[0] == '#':
                continue
            items = line.strip().split()
            genotype = items[3]
            if items[0] == 'duplication':
                Type = 'DUP'
            elif items[0] == 'deletion':
                Type = 'DEL'

            m = re.search('(\w+):(\d+)-(\d+)', items[1])
            chrom = m.group(1)
            start = m.group(2)
            end = m.group(3)
            length = int(end) - int(start)

            Dict[items[3]] = (chrom, start, end, Type, genotype)


    w = open(vcfname, 'w')
    string = '''##fileformat=VCFv4.2
##source=CNVar
##contig=<ID=1,length=249250621>
##contig=<ID=2,length=243199373>
##contig=<ID=3,length=198022430>
##contig=<ID=4,length=191154276>
##contig=<ID=5,length=180915260>
##contig=<ID=6,length=171115067>
##contig=<ID=7,length=159138663>
##contig=<ID=8,length=146364022>
##contig=<ID=9,length=141213431>
##contig=<ID=10,length=135534747>
##contig=<ID=11,length=135006516>
##contig=<ID=12,length=133851895>
##contig=<ID=13,length=115169878>
##contig=<ID=14,length=107349540>
##contig=<ID=15,length=102531392>
##contig=<ID=16,length=90354753>
##contig=<ID=17,length=81195210>
##contig=<ID=18,length=78077248>
##contig=<ID=19,length=59128983>
##contig=<ID=20,length=63025520>
##contig=<ID=21,length=48129895>
##contig=<ID=22,length=51304566>
##contig=<ID=X,length=155270560>
##contig=<ID=Y,length=59373566>
##contig=<ID=MT,length=16569>
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">
##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Type of approach used to detect SV">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of SV">
##INFO=<ID=LEN,Number=1,Type=String,Description="Length of SV">
##INFO=<ID=VARIANT,Number=1,Type=String,Description="Name of Variant">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
'''
    w.write(string)

    ID = 1

    with open(outputfile) as f:
        for line in f:
            items = line.strip().split()
            if items[3] == '1':
                samplename = items[0].split('/')[-1].replace('.base.depth.txt.gz','')
                genotypes = items[1].split('/')[-2].split('_')

                if 'homo' in genotypes:

                    if genotypes[0] != 'wildtype':
                        chrom, start, end, Type, genotype = Dict[genotypes[0]]
                        w.write('%s\t%s\t%s\t-\t<%s>\t.\t.\tSVTYPE=%s;SVMETHOD=CNVar;END=%s;LEN=%s;VARIANT=%s\tGT\t1/1\n' % (chrom, start, ID, Type, Type, end, length,genotypes[0]))
                        ID += 1
                else:
                    if genotypes[0] != 'wildtype':
                        chrom, start, end, Type, genotype = Dict[genotypes[0]]
                        w.write('%s\t%s\t%s\t-\t<%s>\t.\t.\tSVTYPE=%s;SVMETHOD=CNVar;END=%s;LEN=%s;VARIANT=%s\tGT\t0/1\n' % (chrom, start, ID, Type, Type, end, length,genotypes[0]))
                        ID += 1
                    if genotypes[1] != 'wildtype':
                        chrom, start, end, Type, genotype = Dict[genotypes[1]]
                        w.write('%s\t%s\t%s\t-\t<%s>\t.\t.\tSVTYPE=%s;SVMETHOD=CNVar;END=%s;LEN=%s;VARIANT=%s\tGT\t0/1\n' % (chrom, start, ID, Type, Type, end, length,genotypes[1]))
                        ID += 1
    w.close()

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("usage: python %s <CNVar output> <vcf> <mutation_info>" % sys.argv[0])
        exit(1)
    
    outputfile = sys.argv[1]
    vcfname = sys.argv[2]
    mutation_info = sys.argv[3]

    cnvar2vcf( outputfile, vcfname, mutation_info )