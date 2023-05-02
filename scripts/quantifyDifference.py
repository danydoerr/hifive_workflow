#!/usr/bin/env python

from sys import stdout, stderr, exit
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF

import cooler as clr
import numpy as np

from hic import assimilateMatrices, equiSegmentRegions

#QUANTIFIERS = ('pearson', 'own')


#def quantifyDiffs(mtrx1, mtrx2):
#    """ calculate self-made quantifier for difference. We ignore entries that
#    are None """
#
#    assert mtrx1.shape == mtrx2.shape
#
#    comm = 0
#    diff = 0
#    nones = 0
#    for i in xrange(len(mtrx1)):
#        for j in xrange(len(mtrx1)):
#            v1 = mtrx1[i,j]
#            v2 = mtrx2[i,j]
#            if v1 == None or v2 == None:
#                nones += 1
#            if v1 <=0 and v2 <=0:
#                comm += -max(v1, v2)
#            elif v1 >=0 and v2 > 0:
#                comm += min(v1, v2)
#            diff += max(v1, v2) - min(v1, v2)
#    return comm, diff, float(nones)/(mtrx1.shape[0] * mtrx1.shape[1])


def prepareForComparison(mtrx1, mtrx2):
    """ flatten matrices and fill 0s with random noise"""

    assert mtrx1.shape == mtrx2.shape

    m1 = mtrx1.flatten()
    m2 = mtrx2.flatten()
    _isnan = np.logical_or(np.isnan(m1), np.isnan(m2))
    sel = np.logical_and(np.logical_or(m1 != 0, m2 != 0), _isnan == False)
    v1 = np.sum(m1[sel])
    v2 = np.sum(m2[sel])
    c = np.sum(sel == False)
    r1 = np.random.rand(c)
    r2 = np.random.rand(c)

    # this formula makes the assumption that the sum the completely filled
    # matrix m1 is approximately v1*(1+(c/m1.size))/m1.size
    r1 *= (v1*c)/(np.sum(r1)*np.sum(sel))
    r2 *= (v2*c)/(np.sum(r2)*np.sum(sel))
    m1[sel == False] = r1
    m2[sel == False] = r2
    return m1, m2, c

def quantifyPearson(mtrx1, mtrx2):
    """ calculate pearson correlation of two Hi-C matrices"""

    assert mtrx1.shape == mtrx2.shape
    return np.corrcoef(mtrx1, mtrx2)[0,1]


if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ADHF)
#    parser.add_argument('-t', '--type', type=str, choices=QUANTIFIERS, default =
#            QUANTIFIERS[0], help = 'choose type of quantification')
    parser.add_argument('-r', '--raw_counts', action='store_true',
            help='produce matrix with raw counts, rather than balanced counts')
    parser.add_argument('hic_map', type=file, nargs=2,
            help='Hi-C maps in COOLER format')
    args = parser.parse_args()

    c1 = clr.Cooler(args.hic_map[0].name)
    m1 = c1.matrix(balance=not args.raw_counts)[:]
    bins1 = c1.bins()[:]
    coords1 = zip(bins1['chrom'].values, bins1['start'].values,
            bins1['end'].values)

    c2 = clr.Cooler(args.hic_map[1].name)
    m2 = c2.matrix(balance=not args.raw_counts)[:]
    bins2 = c2.bins()[:]
    coords2 = zip(bins2['chrom'].values, bins2['start'].values,
            bins2['end'].values)
    (m1, m2), common = assimilateMatrices(((m1, coords1, coords1), (m2, coords2,
        coords2)), return_common=True)

    mtrx1 = m1[0]
    mtrx2 = m2[0]

    if mtrx1.shape != mtrx2.shape:
        print >> stderr, 'ERROR: the two matrices have different shapes. Exiting.'
        exit(1)

    if 'baited' in bins1.columns and 'baited' in bins2.columns:
        resolution = m1[1][0][2]-m1[1][0][1]
        sel_baited1 = bins1['baited'].values[common[0][0]]
        sel_baited2 = bins2['baited'].values[common[1][0]]
        baited_mtrx1 = mtrx1[np.ix_(sel_baited1, sel_baited1)]
        baited_mtrx2 = mtrx2[np.ix_(sel_baited2, sel_baited2)]

        # only extract upper triangular matrix, because it's symmetric
        upper_tri = np.triu_indices(baited_mtrx1.shape[0], k=1)
        b1, b2, bc = prepareForComparison(baited_mtrx1[upper_tri],
                baited_mtrx2[upper_tri])
        # subtract zeros from lower triangular matrix including the diagonal
        othere_mtrx1 = mtrx1[sel_baited1, :]
        othere_mtrx2 = mtrx2[sel_baited2, :]
        o1, o2, oc = prepareForComparison(othere_mtrx1, othere_mtrx2)
        corr = quantifyPearson(np.concatenate((b1, o1)), np.concatenate((b2,
            o2)))
        outcols = (corr, '', (bc+oc)/float(b1.size+o1.size))
    else:
        upper_tri = np.triu_indices(mtrx1.shape[0], k=1)
        v1, v2, c = prepareForComparison(mtrx1[upper_tri], mtrx2[upper_tri])
        corr = quantifyPearson(v1, v2)
        outcols = (corr, '', c/float(v1.size))


#    if args.type == 'pearson':
#        outcols = quantifyPearson(mtrx1, mtrx2)
#    elif args.type == 'own':
#        outcols = quantifyDiffs(mtrx1, mtrx2)

    # add column for percentage of outliers
    if 'outlier' in bins1.columns and 'outlier' in bins2.columns:
        o1 = bins1['outlier'].values[common[0][0]]
        o2 = bins2['outlier'].values[common[1][0]]
        outcols += (np.logical_or(o1, o2).sum()/float(o1.size))

    print >> stdout, '\t'.join(map(str, outcols))
