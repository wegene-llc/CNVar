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

def safe_open(infile):
    if infile.endswith('.gz'):
        import gzip
        return gzip.open(infile,'rt')
    else:
        return open(infile)


def calculateSamtoolsDepth(depthfile):
    sum_depth = 0
    average_depth = 0
    count = 0.0

    with safe_open(depthfile) as f:
        for line in f:
            items = line.strip().split()
            if items[3] != 'N':
                count += 1
                sum_depth += int(items[2])

    return sum_depth/count



if __name__ == '__main__':
    if(len(sys.argv) <  2):
        print("usage: python %s <depth>" % sys.argv[0])
        exit(0)

    depthfile = sys.argv[1]

    print(calculateSamtoolsDepth(depthfile))


