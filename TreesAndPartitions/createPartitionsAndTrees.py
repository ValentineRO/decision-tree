from sklearn.model_selection import train_test_split
from sklearn import tree
import numpy as np
import random as rd
import matplotlib.pyplot as plt


def readData(namefile):
    with open("../data/"+namefile+".txt") as f :
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

def createPartition(I,namefile):
    N_test = int(np.ceil(I*0.25))
    N_val = N_test
    N_train = I - 2*N_test
    
    index = [i for i in range(I)]
    rd.shuffle(index)
    
    with open(namefile+".txt",'w+') as f:
        f.write(str(N_train)+"\n")
        f.write(str(N_val)+"\n")
        f.write(str(N_test)+"\n")
        
        for i in range(I):
            f.write(str(index[i])+"\n")
    
    return index[:N_train]

def getTreeDepth(tree):
    D = 0
    nodes = [[0,1]]
    while len(nodes)>0:
        node = nodes.pop()
        if tree.children_left[node[0]] == -1:
            D = max(D,node[1])
        else:
            nodes.append([tree.children_left[node[0]],node[1]+1])
            nodes.append([tree.children_right[node[0]],node[1]+1])
    
    return D-1

def writeTree(tree,J,K,namefile):
    D = getTreeDepth(tree)
    
    N = 2**D-1
    L = 2**D
    
    a = np.zeros((J,N))
    b = np.zeros((N,))
    c = -np.ones((N+L,),dtype=int)
    
    nodes = [[0,0]] # the first number is the number for the structure tree, the second is its number in a complete tree
    while len(nodes)>0:
        node = nodes.pop()
        if tree.children_left[node[0]] == -1:
            c[node[1]] = np.argmax(tree.value[node[0]])
        else:
            a[tree.feature[node[0]],node[1]] = 1
            b[node[1]] = tree.threshold[node[0]]
            c[node[1]] = -1
            nodes.append([tree.children_left[node[0]],2*node[1]+1])
            nodes.append([tree.children_right[node[0]],2*node[1]+2])
            
    with open(namefile+".txt",'w+') as f:
        f.write(str(D)+"\n")
        f.write(str(J)+"\n")
        f.write(str(K)+"\n")
        for t in range(N+L):
            if (t<N):
                for j in range(J):
                    f.write(str(a[j,t])+"\n")
                f.write(str(b[t])+"\n")
                f.write("0.0\n") # ca correspond à epsilon
            f.write(str(c[t])+"\n")
            

def generatePartitionsAndTrees(P,DMAX,namefile):
    X,Y,I,J,K = readData(namefile)
    for p in range(P):
        train_indexes = createPartition(I,"./"+namefile+"/partition"+str(p))
        X_train = X[train_indexes,:]
        Y_train = Y[train_indexes]
        
        Nmin = max(int(np.floor(len(Y_train)*0.05)),1)
        
        for d in range(1,DMAX+1):
            for c in range(d,2**d): # c splits gives c+1 leaves 
                clf = tree.DecisionTreeClassifier(max_depth=d,max_leaf_nodes=c+1,min_samples_leaf=Nmin)
                clf = clf.fit(X_train, Y_train)
                writeTree(clf.tree_,J,K,"./"+namefile+"/CART_part"+str(p)+"_D"+str(d)+"_C"+str(c)+"_Nmin")
                
                clf2 = tree.DecisionTreeClassifier(max_depth=d,max_leaf_nodes=c+1)
                clf2 = clf2.fit(X_train, Y_train)
                writeTree(clf2.tree_,J,K,"./"+namefile+"/CART_part"+str(p)+"_D"+str(d)+"_C"+str(c))
    
#generatePartitionsAndTrees(5,4,"iris")
#generatePartitionsAndTrees(5,4,"wine")
#generatePartitionsAndTrees(5,4,"dermatology")
#generatePartitionsAndTrees(5,4,"blood_donation")
#generatePartitionsAndTrees(5,4,"breast_cancer")
generatePartitionsAndTrees(5,4,"exemple_papier")
