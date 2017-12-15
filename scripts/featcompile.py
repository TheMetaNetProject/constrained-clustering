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
import itertools
import scipy
import scipy.io
import math
import sklearn.metrics
import sklearn.cluster
import sys
import warnings
import time

load = 0 #at what stage should we resume the script?
eps = numpy.finfo(float).eps #machine epsilon

basepath = '/u/metanet/clustering/constrained-clustering/'
voc_size = int(sys.argv[2])

suffix = 'BNC2000'

def main():
    n_clusters = int(sys.argv[1])
    print suffix + ' ' + str(n_clusters)
    vocabpath = basepath+'data/vocab'+suffix+'.txt'
    vocab =[word.strip() for word in open(vocabpath).readlines()]
    vocab = vocab[:voc_size]
    matpath = basepath + 'data/F-rels'+str(voc_size)+suffix+'.mat'
    datapath = basepath+'data/all-BNC-EN.txt'
    print str(n_clusters)
    if load<1:
        fmatrix = matbuild(datapath,vocab)
        scipy.io.savemat(matpath,{'fmatrix':fmatrix})
    else:
        a = scipy.io.loadmat(matpath)
        fmatrix = a['fmatrix']
    for i in range(numpy.size(fmatrix,0)):
        fmatrix[i,:] = fmatrix[i,:]/sum(fmatrix[i,:]+eps)
    if load<2:
        simmatrix = simbuild(fmatrix)
        scipy.io.savemat(basepath+'data/F-sims'+str(voc_size)+suffix+'.txt.mat', {'simmatrix': simmatrix})
    else:
        a = scipy.io.loadmat(basepath+'data/F-sims'+str(voc_size)+suffix+'.txt.mat')
        simmatrix = a['simmatrix']
    outputfile = basepath+'data/optlabels-'+str(n_clusters)+'-'+str(voc_size)+suffix+'.mat'
    if load<3:
        (labels,score) = cluster(simmatrix, n_clusters, outputfile)
    else:
        A = scipy.io.loadmat(outputfile)
        labels = A['labels']
    outpath = basepath + 'data/output'+str(n_clusters)+'-'+str(voc_size)+suffix+'.txt'
    outwrite(vocab,labels, outpath)

def cluster(simmatrix, n_clusters, savefile):
    n_iters = 100
    opt_score = -1e10
    for i in range(n_iters):
        print('iter '+str(i)+'\n')
        labels0 = sklearn.cluster.spectral_clustering(simmatrix, n_clusters=n_clusters, n_init = 10)
        score0 = sklearn.metrics.silhouette_score(simmatrix,labels0,"precomputed")
        counts = numpy.bincount(labels0) - 2000/n_clusters
        pows = numpy.power(counts, 2)
        print(str(score0)+'\t'+str(numpy.sum(pows))+'\n')
        score0 = -numpy.sum(numpy.divide(pows, 1e6))
        if score0>opt_score:
            opt_labels = labels0
            opt_score = score0
            scipy.io.savemat(savefile,{'labels':opt_labels, 'score': opt_score})
    return (opt_labels, opt_score)

def outwrite(vocab,labels, outpath):
    outfile = open(outpath,'w')
    labels1 = []
    for i in labels:
        try:
            labels1.append(i[0])
        except:
            labels1.append(i)
    for i in set(labels1):
        cluster = []
        str1 = ''
        for (word,label) in itertools.izip(vocab,labels):
            if label==i:
                cluster.append(word)
                str1 = str1 + word + ' '
        outfile.write("Cluster "+str(i)+":\t" + str1+'\n')

def matbuild(infilepath, vocab):
    """ turn a text file with sparse text coordinates and counts for relations into a numerical matrix"""
    def buildcontextlist(infilepath):
        context = [line.split()[3] for line in open(infilepath).readlines()]
        context = list(set(context))
        return context
    context = buildcontextlist(infilepath)
    outmatrix = numpy.zeros((len(vocab), len(context)))
    for line in open(infilepath).readlines():
        split = line.split()
        try:
            i = vocab.index(split[2])
            freq = int(split[0])
            j = context.index(split[3])
            outmatrix[i,j] = freq
        except ValueError:
            print split[2]
            pass
    return outmatrix

def simbuild(inmat):
    """ build a symmetric similarity matrix """
    numwords = numpy.size(inmat,0)
    simmat = numpy.zeros((numwords,numwords))
    for i in range(numwords):
        a = inmat[i,:]
        f = open(basepath+'/scripts/j.txt', 'w')
        f.write(str(i))
        f.close()
        for j in range(i+1, numwords):
            b = inmat[j,:]
            simmat[i,j] = math.exp(-jsd(a,b))
    simmat = simmat + numpy.transpose(simmat)
    for i in range(numwords):
        simmat[i,i] = 1
    return simmat

def jsd(x,y): #Jensen-shannon divergence
    with numpy.errstate(all='ignore'): 
        x = numpy.array(x)
        y = numpy.array(y)
        d1 = x*numpy.log2(2*x/(x+y))
        d2 = y*numpy.log2(2*y/(x+y))
        d1[numpy.isnan(d1)] = 0
        d2[numpy.isnan(d2)] = 0
        d = 0.5*numpy.sum(d1+d2)    
    return d

if __name__ == '__main__':
    main()
