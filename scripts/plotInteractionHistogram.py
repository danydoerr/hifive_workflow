#!/usr/bin/env python2

from sys import stdout, stderr, exit
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF, \
        FileType
import logging

import cooler as clr
import numpy as np

import os
if not os.environ.get('DISPLAY', None):
    import matplotlib; matplotlib.use('Agg')
import matplotlib.pylab as plt
from matplotlib import cm

from hic import parseCoords 

LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)

QUANTILE = 0.05


if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('-r', '--raw_counts', action='store_true',
            help = 'output matrix with raw counts, rather than balanced counts')
    parser.add_argument('cool_file', type=file,
            help='Hi-C matrix in COOLER format')
    parser.add_argument('-c', '--chromosome_only', type=str,
            help='restrict the analysis to given chromosme')
    parser.add_argument('-l', '--highlight', nargs='+', type=str,
            help='highlight counts from coordinate')
    args = parser.parse_args()


    # setup logging
    ch = logging.StreamHandler(stderr)
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(logging.Formatter('%(levelname)s\t%(asctime)s\t%(message)s'))
    LOG.addHandler(ch)

    if args.highlight and len(args.highlight) > 2:
        LOG.fatal('cannot highlight more than one region, exiting')
        exit(1)

    LOG.info('loading COOLER %s' %args.cool_file.name)
    c = clr.Cooler(args.cool_file.name)

    mtrx = c.matrix(balance = not args.raw_counts)[:]
    if args.chromosome_only:
        start, end = c.extent(args.chromosome_only)
        mtrx = mtrx[start:end+1,start:end+1]

    t1, t2 = np.quantile(mtrx[np.invert(np.isnan(mtrx))], (QUANTILE/2, 1-QUANTILE/2))

    plt.figure()
    plt.hist(mtrx.flatten(), 200, range=(t1, t2), density=True)
    if args.highlight:
        coord1 = coord2 = args.highlight[0]
        if len(args.highlight) == 2:
            coord2 = args.highlight[1]
        int1 = c.extent(coord1)
        int2 = c.extent(coord2)
        x = mtrx[int1[0]:int1[1]+1,int2[0]:int2[1]+1].mean()
        s = np.sqrt(mtrx[int1[0]:int1[1]+1,int2[0]:int2[1]+1].var())
        plt.axvline(x, color='red')

    plt.savefig(stdout, format='pdf')

    LOG.info('DONE')
