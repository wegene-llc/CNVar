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
import random

def cnvSimulator( in_fq, out_fq, chrom, start, end, Type, ins_start,ins_end, outdir ):
    ins_start = int(ins_start)
    ins_end = int(ins_end)

    outfile = open(out_fq,'w')
    string = ''

    w = open('%s/duplication_insert_position.log' % outdir, 'w')

    with open(in_fq) as f:
        line = f.readline().strip()
        if True:
            outfile.write(line + "\n")

            for line in f:
                string = string + line.strip()
            if Type == 'del':
                head = string[0:int(start)-1]
                tail = string[int(end):]
                newstring = head +  tail
                n = 0
                while n <= len(newstring):
                    outfile.write(newstring[n:n+60] + "\n")
                    n += 60
            if Type == 'dup':
                cnv = string[int(start)-1:int(end)]
                while True:
                    rand = random.randint(int(start),int(end))
                    w.write(out_fq + '\n')
                    w.write('insertion position:%s\n' % rand)
                    newstring = string[0:rand] + cnv +string[rand:]
                    n = 0
                    while n <= len(newstring):
                        outfile.write(newstring[n:n+60] + "\n")
                        n += 60
                    break
        else:
            print('please input' + chrom + '.fa')
    
    w.close()
    
    return 0


if __name__ == '__main__':
    if len(sys.argv) < 7:
        print('python cnv_similator.py <in_fq> <out_fq> <chr> <start> <end> <type> <insert start> <insert end> <outdir>')
        exit(1)

    in_fq = sys.argv[1]
    out_fq = sys.argv[2]
    chrom = sys.argv[3]
    start = sys.argv[4]
    end = sys.argv[5]
    Type = sys.argv[6]
    ins_start = sys.argv[7]
    ins_end = sys.argv[8]
    outdir = sys.argv[9]

    cnvSimulator(in_fq, out_fq, chrom, start, end, Type, ins_start,ins_end, outdir)