#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      Administrator
#
# Created:     23/05/2014
# Copyright:   (c) Administrator 2014
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import numpy
import scipy
import scipy.io
import time
import itertools

eps = numpy.finfo(float).eps #machine epsilon

basepath = '/u/metanet/clustering/constrained-clustering/'
suffix = 'BNC2000'

def main():
    vocabpath = basepath+'data/vocab'+suffix+'.txt'
    constraintpath = basepath +'data/constraints'+suffix+'.txt'
    outpath = basepath + 'data/constraints'+suffix+'.mat'
    vocab =[word.strip() for word in open(vocabpath).readlines()]
    constraints = []
    constraints.append([word.split()[0] for word in open(constraintpath).readlines()])
    constraints.append([word.split()[1] for word in open(constraintpath).readlines()])
    row1, row2 = [],[]
    for word1, word2 in itertools.izip(constraints[0], constraints[1]):
        try:
            a = vocab.index(word1)
            b = vocab.index(word2)
            row1.append(a)
            row2.append(b)
            print a + '\t' + b
        except:
            pass
    A = numpy.matrix([row1, row2])
    A = A.transpose()
    scipy.io.savemat(outpath, {'constraints':A})
if __name__ == '__main__':
    main()
