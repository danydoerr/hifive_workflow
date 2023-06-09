#!/usr/bin/env python2

from sys import stdout, stderr, exit
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF, \
        FileType
from itertools import combinations, product
from os.path import basename
import logging

from h5py import File
import numpy as np


LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)

OUTPUT_FORMATS = ('trv', 'bg2')

def writeMtrxTRV(mtrx, colnames, rownames, title, out):
    print >> out, '%s\tRegions\t%s' %(title, '\t'.join(colnames))
    for i in xrange(mtrx.shape[0]):
        out.write('\t'.join((rownames[i], rownames[i])))
        for j in xrange(mtrx.shape[1]):
            out.write('\t')
            if mtrx[i][j] != None:
                out.write(str(mtrx[i][j]))
        out.write('\n')

    i = len(mtrx)
    while i < len(rownames):
        print >> out, '\t'.join((rownames[i], rownames[i]))
        i += 1

def getEndPos(labels, i):
    if len(labels) - 1 > i and labels[i][0] == labels[i+1][0]:
        return labels[i+1][1]-1
    if i > 0 and labels[i][0] == labels[i-1][0]:
        return 2*labels[i][1] - labels[i-1][1]-1
    return labels[i][1]


def writeMtrxBG2(mtrx, labels, out):
    for i in xrange(mtrx.shape[0]-1):
        chr1, start1 = map(str, labels[i])
        end1 = str(getEndPos(labels, i))
        for j in xrange(i, mtrx.shape[1]):
            if mtrx[i][j] != None:
                chr2, start2 = map(str, labels[j])
                end2 = str(getEndPos(labels, j))
                print >> out, '\t'.join((chr1, start1, end1, chr2, start2,
                    end2, str(mtrx[i][j])))


def hifive2mtrx(data, norm_by_expect=True):
    ''' constructs full matrix from hdf5 data '''

    LOG.info('determining full matrix size')
    # identify size of total matrix
    labels = list()
    startpos = [0]
    chrxs = data['chromosomes'][...]
    for chrx in chrxs:
        pos = data['%s.positions' %chrx][...]
        for p in pos[:, 0]:
            labels.append((chrx, p))
        startpos.append(startpos[-1] + pos.shape[0])
    
    LOG.info('full size: %sx%s' %(len(labels), len(labels)))

    # allocate matrix
    mtrx = np.empty((len(labels), len(labels)), dtype=object)
    
    LOG.info('setting up cis-matrices')
    # insert cis contact counts
    for k, chrx in enumerate(chrxs):
        l = startpos[k+1]-startpos[k]
        counts = data['%s.counts' %chrx][...]
        expected = data['%s.expected' %chrx][...]
        c = startpos[k]
        for z, (i, j) in enumerate(combinations(xrange(l), 2)):
            if counts[z] > 0:
                count = counts[z]
                if norm_by_expect:
                    count /= expected[z]
                mtrx[c+i, c+j] = mtrx[c+j, c+i] = count

    LOG.info('setting up trans-matrices')
    # insert trans contact counts
    for (l, chrx), (k, chry) in combinations(enumerate(chrxs), 2):
        cx = startpos[l]
        cy = startpos[k]
        sx = startpos[l+1]-cx
        sy = startpos[k+1]-cy
        counts = data['%s_by_%s.counts' %(chrx, chry)][...]
        expected = data['%s_by_%s.expected' %(chrx, chry)][...]
        for i, j in product(xrange(sx), xrange(sy)):
            if counts[i,j] > 0:
                count = counts[i,j]
                if norm_by_expect:
                    count /= expected[i,j]
                mtrx[cx+i,cy+j] = mtrx[cy+j,cx+i] = count
        
    return mtrx, labels


if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('hdf5_file', type=file, 
            help='Hi-C matrix in HDF5 format generated by HiFive')
    parser.add_argument('-b', '--output_bed', type = FileType('w'),
            help='output segmentation of the genome in BED (version 2) format')
    parser.add_argument('-n', '--normalize_by_expected', action='store_true',
            help='normalize (divide) counts by their expected value')
    parser.add_argument('-f', '--out_format', type=str, choices=OUTPUT_FORMATS,
            default=OUTPUT_FORMATS[0], help='choose output format of matrix')
    args = parser.parse_args()

    # setup logging
    ch = logging.StreamHandler(stderr)
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(logging.Formatter('%(levelname)s\t%(asctime)s\t%(message)s'))
    LOG.addHandler(ch)

    LOG.info('reading HDF5 file %s' %args.hdf5_file.name)
    data = File(args.hdf5_file.name, 'r')
    LOG.info('file contains data from %s chromosomes' %(data['chromosomes'].len()))

    LOG.info('transforming into single matrix')
    mtrx, labels = hifive2mtrx(data, args.normalize_by_expected)
    if args.out_format == 'trv':
        LOG.info('writing output in TRV format')
        names = map(lambda x: '%s-%s' %x, labels)
        writeMtrxTRV(mtrx, names, names, '%s:count%s' %(basename( \
                args.hdf5_file.name), args.normalize_by_expected and \
                '/expected' or ''), stdout)
    elif args.out_format == 'bg2':
        LOG.info('writing output in BG2 format')
        writeMtrxBG2(mtrx, labels, stdout)
    if args.output_bed:
        for i, (chrx, startx) in enumerate(labels):
            end = 0
            if i + 1 < len(labels) and labels[i+1][0] == chrx:
                end = labels[i+1][1]
            else:
                end = 2*labels[i][1]-labels[i-1][1]
            print >> args.output_bed, '\t'.join((chrx, str(labels[i][1]),
                str(end)))
        args.output_bed.close()

    LOG.info('DONE')
