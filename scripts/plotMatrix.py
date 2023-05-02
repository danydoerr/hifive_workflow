#!/usr/bin/env python

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

LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)


if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('-c', '--raw_counts', action='store_true',
            help = 'output matrix with raw counts, rather than balanced counts')
    parser.add_argument('cool_file', type=file,
            help='Hi-C matrix in COOLER format')
    args = parser.parse_args()

    # setup logging
    ch = logging.StreamHandler(stderr)
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(logging.Formatter('%(levelname)s\t%(asctime)s\t%(message)s'))
    LOG.addHandler(ch)

    LOG.info('loading COOLER %s' %args.cool_file.name)
    c = clr.Cooler(args.cool_file.name)

    mtrx = c.matrix(balance = not args.raw_counts)[:]

    #mx = np.quantile(mtrx[np.isfinite(mtrx)], 0.99)
    mx = 0.00045
    mtrx[mtrx > mx] = mx
    m = mtrx * len(cm.viridis.colors)/float(mx)
    axMtrx = plt.axes([0.95, 0, 0.95, 0.95], frameon=False)
    axMtrx.imshow(m, cmap=cm.viridis, interpolation='none')
    axMtrx.axis('off')
    axMtrx.axes.get_xaxis().set_visible(False)
    axMtrx.axes.get_yaxis().set_visible(False)
    axColor = plt.axes([0.95, 0, 0.1, 0.95], frameon=False)
    axColor.imshow(plt.transpose(( plt.linspace(1, 0, 256),)),
            plt.get_cmap('viridis'),
        aspect='auto', extent=[0,1,0, mx] )
    axColor.set_xticks([])

    plt.savefig(stdout, format='pdf', bbox_inches="tight", pad_inches=0.1,
            dpi=300)

    LOG.info('DONE')
