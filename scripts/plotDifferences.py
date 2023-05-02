#!/usr/bin/env python

from sys import stdout, stderr, exit
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF
from os.path import basename
import re

import numpy as np
import os
if not os.environ.get('DISPLAY', None):
    import matplotlib; matplotlib.use('Agg')
import matplotlib.pylab as plt


if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('-p', '--pearson', action='store_true',
            help='indicate that input is pearson correlation, not own measure')
    parser.add_argument('difference_stats', type=file,
            help='difference stats of same probes, but different genomic regions')
    args = parser.parse_args()

    name = basename(args.difference_stats.name).rsplit('.', 1)[0]
    stats = plt.loadtxt(args.difference_stats, ndmin=2)

    min_ints = sorted(set(stats[:, 0]))
    resolution = sorted(set(stats[:, 1]))
    mxcols = 6

    nrows = 1+(len(min_ints)-1)/mxcols
    ncols = min(mxcols, len(min_ints))
    fig = plt.figure(figsize=(ncols * 4, nrows * 4 * 2), dpi=300)

    spi = 0
    superres = None
    if args.pearson:
        sim = lambda x: x[:, 1]
    else:
        sim = lambda x: x[:, 1]/(x[:, 1]+x[:, 2])
    for m in min_ints:
        spi+=1
        plt.subplot(nrows * 2, ncols, spi)
        for r in resolution:
            c = stats[(stats[:, 0] == m) & (stats[:, 1] == r)][:, 2:]
            p = plt.plot(c[:, 0], sim(c), label='%skb' %int(r/1000))
            #p = plt.plot(c[:, 0], sim(c), label='%skb' %int(r/1000))
        #plt.xlim([resolution[0], resolution[-1]])
        #plt.xlabel('super-resolution')
        if not ((spi-1) % ncols):
            plt.ylabel('similarity')
            if args.pearson:
                plt.ylabel('Pearson correlation')
        _sims = sim(stats[:, 2:])
        _sims_not_nan = _sims[np.isnan(_sims) == False]
        plt.ylim([min(_sims_not_nan)*0.9, max(_sims_not_nan)*1.1])
#        plt.ylim([min(sim(stats[:, 2:]))*0.9, max(sim(stats[:, 2:]))*1.1])
        #plt.axhline(y=0.9, color='r', linestyle='--')
        #plt.ylim([0, 1])
        plt.legend(loc='lower right')
        plt.title('min interactions %s' %m, fontsize=10)

    spi = nrows * ncols
    for m in min_ints:
        spi+=1
        plt.subplot(nrows * 2, ncols, spi)
        for r in resolution:
            c = stats[(stats[:, 0] == m) & (stats[:, 1] == r)][:, 2:]
            if c.shape[1] >= 4:
                p = plt.plot(c[:, 0], c[:, 2]*100, label='%skb empty cells' %int(r/1000))
                p = plt.plot(c[:, 0], c[:, 3]*100, '--', color=p[0].get_color(),
                        label='%skb outlier ' %int(r/1000))
            else:
                p = plt.plot(c[:, 0], c[:, 2]*100, label='%skb' %int(r/1000))
        #plt.xlim([resolution[0], resolution[-1]])
        plt.xlabel('super-resolution')
        if not ((spi-1) % ncols):
            plt.ylabel('percentage')
        #plt.axhline(y=0.9, color='r', linestyle='--')
        data = stats[:, 4]
        plt.ylim([0,  max(data[np.isnan(data) == False].flatten())*110])
        if stats.shape[1] < 7:
            plt.title('empty cells')

        plt.legend(loc='upper right')
        #plt.title('min interactions %s' %m, fontsize=10)

    fig.suptitle('difference between %s' %name)

#    for i, r in enumerate(resolution):
#        spi += 1
#        plt.subplot(nrows, ncols, spi)
#        c = stats[(stats[:, 0] == m) & (stats[:, 1] == r)][:, 2:]
#        h = plt.boxplot([[f(stats[(stats[:, 0] == r) & (stats[:, 1] == s)]) for _,
#            stats in data] for s in superres], labels=map(lambda
#                x: str(int(x)), superres))
#        for el in ['boxes', 'whiskers', 'fliers', 'means', 'medians', 'caps']:
#            plt.setp(h[el], color=colors[i])
#
#        plt.axhline(y=0.9, color='r', linestyle='--')
#        plt.ylim([0, 1])
#        plt.ylabel('similarity')
#        plt.xlabel('super-resolution')
#        plt.title('%skb' %int(r/1000), fontsize=6)

    #plt.show()
    plt.subplots_adjust(top = 0.9, bottom=0.1, hspace=0.25, wspace=0.3)
    plt.savefig(stdout, format='pdf', bbox_inches="tight", pad_inches=0.1)

