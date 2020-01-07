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

def trans_seq(seq):
    def trans_base(base):
        if base == 'T':
            return 'A'
        elif base == 'A':
            return 'T'
        elif base == 'C':
            return 'G'
        elif base == 'G':
            return 'C'
        elif base == 'N':
            return 'N'

    def trans_loc(seq):
        seq2 = ""
        for i in range(len(seq)):
            seq2 += trans_base(seq[-(i+1)])
        return seq2

    return trans_loc(seq)


def kMers(infile, outfile1, outfile2, k, insert_size):
    # Read infile
    header = 'HEADER'
    string = ''
    k = int(k)
    insert_size = int(insert_size) - k

    with open(infile) as f:
        for line in f:
            if line[0] == '>':
                header = line.strip().replace(' ',':')
            else:
                string += line.strip()

    fq1 = open(outfile1,'w')
    fq2 = open(outfile2,'w')

    for pos in range(0, len(string)-k+1):
        fq1.write('%s:%s\n' %(header,pos))
        fq1.write(string[pos:pos+k] + '\n')
        fq2.write('%s:%s\n' %(header,pos))
        fq2.write(trans_seq(string[pos+insert_size:pos+k+insert_size]) + '\n')

    return 0

if __name__ == '__main__':
    if len(sys.argv) < 5:
        print("usage: python k-mers.py <in.fasta> <out.fasta1> <out.fasta2> <k value> <insert size>")
        exit(1)

    infile = sys.argv[1]
    outfile1 = sys.argv[2]
    outfile2 = sys.argv[3]
    k = sys.argv[4]
    insert_size = sys.argv[5]

