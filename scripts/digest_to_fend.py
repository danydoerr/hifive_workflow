#!/usr/bin/env python

from sys import stdout, stderr, exit
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF, \
        FileType
import csv


if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('digest_file', type=file, 
            help='tab-separated list of restriction digest fragments ')
    args = parser.parse_args()
    
    out = stdout

    header = 2
    c = 1
    print >> out, '\t'.join(('fend', 'frag', 'chr', 'coord', 'valid',
        'frag_len'))
    for line in csv.reader(args.digest_file, delimiter='\t'):
        if header:
            header -=1
            continue

        flen = str(int(line[2])-int(line[1])+1)
        print >> out, '\t'.join((str((c-1)*2+1), str(c), line[0], line[1], '1', flen))
        print >> out, '\t'.join((str(c*2), str(c), line[0], line[2], '1', flen))

        c += 1 
