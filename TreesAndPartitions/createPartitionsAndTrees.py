from sklearn import tree
import numpy as np
import random as rd
import pandas as pd
import matplotlib.pyplot as plt

nbOfPart = 5

dtNamesWithTrain = ["car_evaluation","IndiansDiabetes","Ionosphere","iris","messidor","monk1","monk2","monk3","seismic_bumps",
                    "spambase","Statlog_satellite","tic-tac-toe","wine"]
dtNamesWithoutTrain = ["balance-scale","bank_conv","banknote","biodeg"]
dtNamesWithoutPart = ["breast_cancer","blood_donation","dermatology","digits","ecoli","german","haberman","seeds"]

def establishPartitions(name, withTrain, prec = 0, sep = ";"):
    if name == "Ionosphere":
        dataset = pd.read_csv("./data/datasets_SiccoVerwer/"+name+".csv", sep = ',', dtype="float")
    else:
        dataset = pd.read_csv("./data/datasets_SiccoVerwer/"+name+".csv", sep = sep, dtype="float")
    
    I = dataset.shape[0]
    Jp = dataset.shape[1]
    
    cols = dataset.columns
    
    for p in range(5):
        if withTrain:
            partTrain = pd.read_csv("./data/datasets_SiccoVerwer/"+name+".csv.train"+str(p)+".csv", sep = sep, dtype="float")
            partTest = pd.read_csv("./data/datasets_SiccoVerwer/"+name+".csv.test"+str(p)+".csv", sep = sep, dtype="float")
            
            I_train, I_val, I_test = [], [], []
            
            for i in range(I):
                df = partTrain[partTrain[cols[0]] == dataset.iloc[i][0]] 
                for j in range(1,Jp):
                    if df.shape[0] == 0:
                        break
                    else:
                        if prec !=0:
                            df = df[(df[cols[j]] > dataset.iloc[i][j]-prec)&(df[cols[j]] < dataset.iloc[i][j] + prec)]
                        else:
                            df = df[df[cols[j]] == dataset.iloc[i][j]]
                if df.shape[0] > 0:
                    partTrain = partTrain.drop(index=df.iloc[0].name)
                    I_train.append(i)
                else:
                    df = partTest[partTest[cols[0]] == dataset.iloc[i][0]] 
                    for j in range(1,Jp):
                        if df.shape[0] == 0:
                            break
                        else:
                            if prec !=0:
                                df = df[(df[cols[j]] > dataset.iloc[i][j]-prec)&(df[cols[j]] < dataset.iloc[i][j] + prec)]
                            else:
                                df = df[df[cols[j]] == dataset.iloc[i][j]]
                    if df.shape[0] > 0:
                        partTest = partTest.drop(index=df.iloc[0].name)
                        I_test.append(i)
                    else:
                        I_val.append(i)
            
            with open("./"+name+"/"+"partition"+str(p)+".txt",'w+') as f:
                f.write(str(len(I_train))+"\n")
                f.write(str(len(I_val))+"\n")
                f.write(str(len(I_test))+"\n")
                
                for i in range(len(I_train)):
                    f.write(str(I_train[i])+"\n")
        
                for i in range(len(I_val)):
                    f.write(str(I_val[i])+"\n")
        
                for i in range(len(I_test)):
                    f.write(str(I_test[i])+"\n")
        else:
            partTest = pd.read_csv("./data/datasets_SiccoVerwer/"+name+".csv.test"+str(p)+".csv", sep = sep, dtype="float")
        
            I_train, I_test = [], []
            
            for i in range(I):
                df = partTest[partTest[cols[0]] == dataset.iloc[i][0]] 
                for j in range(1,Jp):
                    if df.shape[0] == 0:
                        break
                    else:
                        if prec !=0:
                            df = df[(df[cols[j]] > dataset.iloc[i][j]-prec)&(df[cols[j]] < dataset.iloc[i][j] + prec)]
                        else:
                            df = df[df[cols[j]] == dataset.iloc[i][j]]
                if df.shape[0] > 0:
                    partTest = partTest.drop(index=df.iloc[0].name)
                    I_test.append(i)
                else:
                    I_train.append(i)
                    
            IVal = int(np.floor(len(I_train)/3))
            ITr = len(I_train) - IVal
            rd.shuffle(I_train)
            
            with open("./"+name+"/"+"partition"+str(p)+".txt",'w+') as f:
                f.write(str(ITr)+"\n")
                f.write(str(IVal)+"\n")
                f.write(str(len(I_test))+"\n")
                
                for i in range(ITr):
                    f.write(str(I_train[i])+"\n")
        
                for i in range(ITr,len(I_train)):
                    f.write(str(I_train[i])+"\n")
        
                for i in range(len(I_test)):
                    f.write(str(I_test[i])+"\n")
        
    createDatasetFile(name)
    
def createDerm():
    with open("./data/UCI/dermatology.data","r") as f:
        X = []
#        Y = []
        line = f.readline()
        while(line != ''):
            X.append(line.split(','))
            if '?' in X[-1]:
                print(X.pop())
#            Y.append(X[-1].pop())
            line = f.readline()
    
    I = len(X)
    Jp = len(X[0])
    J = Jp - 1
    
    cols = ["Col"+str(i) for i in range(Jp)]
    
    dataset = pd.DataFrame(np.array(X), dtype=float, columns = cols)
    
    att_min, att_max = [], []
    for j in range(J):
        att_min.append(np.min(dataset[cols[j]].tolist()))
        att_max.append(np.max(dataset[cols[j]].tolist()))
    
    lab = sorted(dataset[cols[-1]].unique())
    
    with open("./data/dermatology.txt","w+") as f:
        f.write(str(I)+"\n")
        f.write(str(J)+"\n")
        f.write(str(len(lab))+"\n")
        for i in range(I):
            for j in range(J):
                xij = (dataset.iloc[i][j] - att_min[j])/(att_max[j] - att_min[j])
                f.write(str(xij)+"\n")
        for i in range(I):
            yi = lab.index(dataset.iloc[i][-1])
            f.write(str(yi)+"\n")
            
def createBloodT():
    with open("./data/UCI/transfusion.data","r") as f:
        X = []
#        Y = []
        line = f.readline()
        while(line != ''):
            X.append(line.split(','))
#            Y.append(X[-1].pop())
            line = f.readline()
    
    I = len(X)
    Jp = len(X[0])
    J = Jp - 1
    
    cols = ["Col"+str(i) for i in range(Jp)]
    
    dataset = pd.DataFrame(np.array(X), dtype=float, columns = cols)
    
    att_min, att_max = [], []
    for j in range(J):
        att_min.append(np.min(dataset[cols[j]].tolist()))
        att_max.append(np.max(dataset[cols[j]].tolist()))
    
    lab = sorted(dataset[cols[-1]].unique())
    
    with open("./data/blood-transfusion.txt","w+") as f:
        f.write(str(I)+"\n")
        f.write(str(J)+"\n")
        f.write(str(len(lab))+"\n")
        for i in range(I):
            for j in range(J):
                xij = (dataset.iloc[i][j] - att_min[j])/(att_max[j] - att_min[j])
                f.write(str(xij)+"\n")
        for i in range(I):
            yi = lab.index(dataset.iloc[i][-1])
            f.write(str(yi)+"\n")

def createBreastC():
    with open("./data/UCI/breast-cancer-wisconsin.data","r") as f:
        X = []
#        Y = []
        line = f.readline()
        while(line != ''):
            X.append(line.split(','))
            X[-1] = X[-1][1:] # the first "itemé is ID number
            if '?' in X[-1]:
                print(X.pop())
#            Y.append(X[-1].pop())
            line = f.readline()
    
    I = len(X)
    Jp = len(X[0])
    J = Jp - 1
    
    cols = ["Col"+str(i) for i in range(Jp)]
    
    dataset = pd.DataFrame(np.array(X), dtype=float, columns = cols)
    
    att_min, att_max = [], []
    for j in range(J):
        att_min.append(np.min(dataset[cols[j]].tolist()))
        att_max.append(np.max(dataset[cols[j]].tolist()))
    
    lab = sorted(dataset[cols[-1]].unique())
    
    with open("./data/breast-cancer.txt","w+") as f:
        f.write(str(I)+"\n")
        f.write(str(J)+"\n")
        f.write(str(len(lab))+"\n")
        for i in range(I):
            for j in range(J):
                xij = (dataset.iloc[i][j] - att_min[j])/(att_max[j] - att_min[j])
                f.write(str(xij)+"\n")
        for i in range(I):
            yi = lab.index(dataset.iloc[i][-1])
            f.write(str(yi)+"\n")
        

def createDatasetFile(name, sep = ";"):
    if name == "Ionosphere":
        dataset = pd.read_csv("./data/datasets_SiccoVerwer/"+name+".csv", sep = ',', dtype="float")
        del dataset["Feat1"]
    elif name == "seismic_bumps":
        dataset = pd.read_csv("./data/datasets_SiccoVerwer/"+name+".csv", sep = sep, dtype="float")
        del dataset["Feat13"]
        del dataset["Feat14"]
        del dataset["Feat15"]
    else:
        dataset = pd.read_csv("./data/datasets_SiccoVerwer/"+name+".csv", sep = sep, dtype="float")
        
    I = dataset.shape[0]
    Jp = dataset.shape[1]
    J = Jp - 1
    
    cols = dataset.columns
    
    att_min, att_max = [], []
    for j in range(J):
        att_min.append(np.min(dataset[cols[j]].tolist()))
        att_max.append(np.max(dataset[cols[j]].tolist()))
    
    lab = sorted(dataset[cols[-1]].unique())
    
    with open("./data/"+name+".txt","w+") as f:
        f.write(str(I)+"\n")
        f.write(str(J)+"\n")
        f.write(str(len(lab))+"\n")
        for i in range(I):
            for j in range(J):
                xij = (dataset.iloc[i][j] - att_min[j])/(att_max[j] - att_min[j])
                f.write(str(xij)+"\n")
        for i in range(I):
            yi = lab.index(dataset.iloc[i][-1])
            f.write(str(yi)+"\n")

def createPartition(namefile,I):
    for p in range(5):
        N_test = int(np.ceil(I*0.25))
        N_train = int(np.ceil(I*0.5))
        N_val = I - N_test - N_train
        
        if p==0:
            print(N_train)
            print(N_val)
            print(N_test)
        
        index = [i for i in range(I)]
        rd.shuffle(index)
        
        with open("./"+namefile+"/partition"+str(p)+".txt",'w+') as f:
            f.write(str(N_train)+"\n")
            f.write(str(N_val)+"\n")
            f.write(str(N_test)+"\n")
            
            for i in range(I):
                f.write(str(index[i])+"\n")

def readPartition(namefile):
    I_train = []
    I_val = []
    I_test = []
    
    with open(namefile+".txt") as f :
        N_train = int(f.readline())
        N_val = int(f.readline())
        N_test = int(f.readline())
        
        for i in range(N_train):
            I_train.append(int(f.readline()))
            
        for i in range(N_val):
            I_val.append(int(f.readline()))
            
        for i in range(N_test):
            I_test.append(int(f.readline()))
            
    return I_train, I_val, I_test

def viewPartitions(namefile, sep =";"):
    if namefile == "Ionosphere":
        dataset = pd.read_csv("./data/datasets_SiccoVerwer/"+namefile+".csv", sep = ',', dtype="float")
    else:
        dataset = pd.read_csv("./data/datasets_SiccoVerwer/"+namefile+".csv", sep = sep, dtype="float")
    
    I = dataset.shape[0]
    
    A = np.zeros((nbOfPart*3,I))
    
    for p in range(nbOfPart):
        Itr,Iv,Itst = readPartition("./"+namefile+"/partition"+str(p))
        
        for i in range(len(Itr)):
            A[p*3,Itr[i]] = 1
    
        for i in range(len(Iv)):
            A[p*3+1,Iv[i]] = 1
    
        for i in range(len(Itst)):
            A[p*3+2,Itst[i]] = 1

        print("Partition ",p)
        print("Taille train: ",len(Itr))
        print("Taille validation: ",len(Iv))
        print("Taille test: ",len(Itst))
        
    plt.close()
    plt.imshow(A,cmap="binary")

def readData(namefile):
    with open("./data/"+namefile+".txt") as f :
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

def createTrees(namefile, DMAX):
    X,Y,I,J,K = readData(namefile)
    for p in range(nbOfPart):
        Itr, Ival, Itst = readPartition("./"+namefile+"/partition"+str(p))
        
        X_train = X[Itr,:]
        Y_train = Y[Itr]
        
        Nmin = max(int(np.floor(len(Y_train)*0.05)),1)
        
        for d in range(1,DMAX+1):
            for c in range(d,2**d): # c splits gives c+1 leaves 
                clf = tree.DecisionTreeClassifier(max_depth=d,max_leaf_nodes=c+1,min_samples_leaf=Nmin)
                clf = clf.fit(X_train, Y_train)
                writeTree(clf.tree_,J,K,"./"+namefile+"/CART_part"+str(p)+"_D"+str(d)+"_C"+str(c)+"_Nmin")
                
                clf2 = tree.DecisionTreeClassifier(max_depth=d,max_leaf_nodes=c+1)
                clf2 = clf2.fit(X_train, Y_train)
                writeTree(clf2.tree_,J,K,"./"+namefile+"/CART_part"+str(p)+"_D"+str(d)+"_C"+str(c))

#datasets = ["car_evaluation","IndiansDiabetes","Ionosphere","iris","messidor","monk1","monk2","monk3","seismic_bumps",
#                    "spambase","Statlog_satellite","tic-tac-toe","wine","balance-scale","bank_conv","biodeg",
#                    "breast_cancer","blood_donation","dermatology","digits","ecoli","german","haberman","seeds"]

datasets = ["chess","HTRU_2","magic04","winequality-red"]

for d in range(len(datasets)):
    print(datasets[d])
    createTrees(datasets[d],4)