#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from sklearn.cluster import AgglomerativeClustering

def readDt():
    with open("./dataset.txt") as f :
        I = int(f.readline())
        J = int(f.readline())
        K = int(f.readline())
        
        X = np.zeros((I,J))
        Y = np.zeros((I,1))
        
        for i in range(I):
            for j in range(J):
                X[i,j] = float(f.readline())
        
        for i in range(I):
            Y[i] = int(f.readline())
    return X,Y,I,J,K

def readParameters():
    with open("./parameters.txt") as f:
        gam = float(f.readline())
    return gam

def transformDt(X,Y,J,K): # dataset of J charact: add K-1 charact to do hierarchical clustering
    I = len(Y)
    
    Xp = np.zeros((I,K-1))
    
    for i in range(I):
        for k in range(K-1):
            Xp[i,k] = J*(Y[i] == k)
    
    Xc = np.concatenate((X,Xp),axis=1)
    return Xc

def writeCl(C):
    with open("./HCL.txt",'w+') as f:
        for c in range(len(C)):
            f.write(str(len(C[c]))+"\n")
            for i in range(len(C[c])):
                f.write(str(C[c][i])+"\n")
        f.write(str(-1))
        
#def showRep(C):
#    nbOfSize = {len(C[0]):1}
#    for c in range(1,len(C)):
#        if len(C[c]) in nbOfSize:
#            nbOfSize[len(C[c])] += 1
#        else:
#            nbOfSize[len(C[c])] = 1
#    
#    for i, (k, v) in enumerate(nbOfSize.items()):
#        print(v, " clusters of size ",k)
        
X,Y,I,J,K = readDt()
Xc = transformDt(X,Y,J,K)

gam = readParameters()
nbCl = max(K,int(len(Y)*gam))

hierarchical_cluster = AgglomerativeClustering(n_clusters=nbCl, affinity='euclidean', linkage='single')
labels = hierarchical_cluster.fit_predict(Xc)

C = []
for c in range(nbCl):
    C.append([])
    for i in range(len(Y)):
        if labels[i] == c:
            C[-1].append(i)
                
writeCl(C)
#showRep(C)
