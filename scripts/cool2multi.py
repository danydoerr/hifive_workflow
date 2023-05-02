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
    parser.add_argument('-c', '--raw_counts', action='store_true',
            help = 'output matrix with raw counts, rather than balanced counts')
    parser.add_argument('-s', '--chromosome_sizes', type=FileType('w'),
            help = 'write chromosome sizes to specified file')
    parser.add_argument('-o', '--out_file', type=FileType('w'),
            default='out.multi.cool',
            help='Output file in COOLER format')
    parser.add_argument('cool_file', type=str, nargs='+',
            help='Hi-C matrix in COOLER format')
    args = parser.parse_args()

    # setup logging
    ch = logging.StreamHandler(stderr)
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(logging.Formatter('%(levelname)s\t%(asctime)s\t%(message)s'))
    LOG.addHandler(ch)

    for f in args.cool_file:
        LOG.info('loading COOLER %s' %f)
        c = clr.Cooler(f)
        pixels = c.pixels()[:]
        bins = c.bins()[:]
        fst_row = bins.head(1)
        resolution = int(fst_row['end']) - int(fst_row['start'])
        uri = '%s::/resolutions/%s' %(args.out_file.name, resolution)

        LOG.info('writing COOLER to %s' %uri)
        clr.create_cooler(uri, bins=bins, pixels=pixels, ordered=True,
                assembly='tair10', mode='a', dtypes = dict((('bin1_id', int),
                ('bin2_id', int), ('count', float))))

#    if args.chromosome_sizes:
#        LOG.info('writing chromosome sizes')
#        for name, size in zip(c.chromnames, c.chromsizes):
#            print >> args.chromosome_sizes, '\t'.join((name, str(size)))
#        args.chromosome_sizes.close()

    LOG.info('DONE')
