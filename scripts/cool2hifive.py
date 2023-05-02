#!/usr/bin/env python

from sys import stdout, stderr, exit
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF, \
        FileType
from itertools import repeat, combinations, chain, product, compress
from os.path import basename
import logging
import csv

import cooler as clr
import numpy as np

LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)


if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('-o', '--out_file', type=FileType('w'),
            default='<input_file>.mat',
            help='Output file in COOLER format')
    parser.add_argument('cool_file', type=str,
            help='Hi-C matrix in COOLER format')
    args = parser.parse_args()

    output_name = args.out_file
    if parser.get_default('out_file') == args.out_file:
        output_name = args.cool_file.rsplit('.', 1)[0]

    # setup logging
    ch = logging.StreamHandler(stderr)
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(logging.Formatter('%(levelname)s\t%(asctime)s\t%(message)s'))
    LOG.addHandler(ch)

    LOG.info('loading COOLER %s' %args.cool_file)
    c = clr.Cooler(args.cool_file)
    pixels = c.pixels()[:]
    bins = c.bins()[:]
    fst_row = bins.head(1)
    resolution = int(fst_row['end']) - int(fst_row['start'])

    fend_out = open('%s.txt' %output_name, 'w')
    LOG.info('writing FEND raw data to %s' %fend_out.name)
    print >> fend_out, '\t'.join(('fend', 'frag', 'chr', 'coord', 'valid',
        'frag_len'))
    for i, line in enumerate(bins.values):
        print >> fend_out, '\t'.join(map(str, (2*i+1, i+1, line[0], line[1]+1, 1,
            resolution)))
        print >> fend_out, '\t'.join(map(str, (2*i+2, i+1, line[0], line[2], 1,
            resolution)))
    fend_out.close()

    mat_out = open('%s.txt' %output_name, 'w')
    LOG.info('writing MAT data to %s' %mat_out.name)
    for i, j, c in pixels.values:
        print >> mat_out, '\t'.join((i+1, j+1, c))
    mat_out.close()
    LOG.info('DONE')
