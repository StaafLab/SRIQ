#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 10:19:05 2020

@author: jacob
"""
from pathlib import Path
import os
from itertools import combinations
from collections import Counter
from collections import defaultdict
import pandas as pd
import numpy as np
from os import listdir
from os.path import isfile, join
import seaborn as sns
import matplotlib.pyplot as plt
np.random.seed(10)

def runSRIQ(data, studyName = 'SRIQ', cutOff = None, permutations = 10000, iterations = 10, minBagSize=1200, minClusterSize = 0):
    #data: string with the directory to the csv file which you want to run
    #cutOff: list of the cutoffs
    output = ''
    resources = '../software/VRLA/resources/test.properties'
    if cutOff is None: cutOff = [0.9,0.89,0.88,0.87,0.86,0.85,0.84,0.83,0.82,0.81,0.80,0.79,0.78,0.77,0.76,0.75,0.74,0.73,0.72,0.71,0.7,0.69,0.68,0.67,0.66,0.65,0.64,0.63,0.62,0.61,0.6,0.59,0.58,0.57,0.56,0.55,0.54,0.53,0.52,0.51,0.5,0.49,0.48,0.47,0.46,0.45,0.44,0.43,0.42,0.41,0.4,0.39,0.38,0.37,0.36,0.35,0.34]
    with open(resources) as file:
        for line in file.readlines():
            if 'studyName' in line:
                output += f'studyName={studyName}\n'
            elif 'inFileName' in line:
                output += f'inFileName={data}\n'
            elif 'distCutOff' in line:
                if isinstance(cutOff, list): 
                    output += 'distCutOff={}\n'.format(", ".join([str(x) for x in cutOff]))
                else:
                    output += f'distCutOff={str(cutOff)}\n'
            elif 'permutations' in line:
                output += f'permutations={permutations}\n'
            elif 'minClusterSize' in line:
                output += f'minClusterSize={minClusterSize}\n'
            elif 'minBagSize' in line:
                output += f'minBagSize={minBagSize}\n'
            elif 'iterations' in line:
                output += f'iterations={iterations}\n'
            else:
                output += line
    with open(resources, 'w') as file:
        file.write(output)
    print('running...')
    bashCommand = "cd ../software/VRLA && java -jar -Xms12g VRLA.jar"
    os.system(bashCommand)
    print('Done!')
    
def splitRunSRIQ(data, frac, studyName, minClusterSize = 0):
    data = f'data/expressionData/{data}'
    df = pd.read_csv(data, sep = '\t')
    data = 'data/expressionData/filtered(21k).txt'
    columns = df.columns.to_list()
    remove_n = int(df.shape[1]*(frac))
    keep = np.random.choice(df.columns, remove_n, replace=False)
    keep = np.insert(columns, 0, 'Gene')
    columns = [ele for ele in columns if ele in keep]
    df = df[columns]
    df = df.drop_duplicates(['Gene'], keep ='first')
    df = df.loc[:,~df.columns.duplicated()]
    print(df.shape)
    df.to_csv(f'data/expressionData/splitData.txt', sep = '\t', index = False)
    runSRIQ('splitData', studyName = studyName, minClusterSize = minClusterSize)



def robustnessAlgorithm(data = 'test_mc_log2var(80)', cutoff = None, permutations = (10000, ), iterations = (10,), minBagSizes = (200,400, 600), minClusterSizes = (0,10, 20), cluster = 4):
    collections = defaultdict(list)
    size = 0
    for permutation in permutations:
        for iteration in iterations:
            for minBagSize in minBagSizes:
                for minClusterSize in minClusterSizes:
                    size += 1
                    print(f'Iteration {size}')
                    runSRIQ(data=data, permutations = permutation, iterations = iteration, minBagSize=minBagSize, minClusterSize = minClusterSize)
                    path = f'/home/jacob/examensarbete/software/output/VRLA_test_{permutation}itr_{minBagSize}var_{iteration}r/{permutation}/QC_Spiral(false)/dist(0.42)/{cluster}'
                    for counter, filename in enumerate(sorted(os.listdir(path))):
                        with open(os.path.join(path, filename)) as file:
                            for x in combinations([x.replace('\n', '') for x in file.readlines()[1:]], 1):
                                collections[counter+1].append(x)
    
    return [(sum(Counter(collections[i+1]).values())/len(Counter(collections[i+1]))/size) for i in range(cluster)]

def parameterRobustness(data = 'test_mc_log2var(80)', cutoff = None, permutations = (10000, ), iterations = (10,), minBagSizes = (200,400, 600), minClusterSizes = (0,)):
    collections = defaultdict(list)
    size = 0
    for minBagSize in minBagSizes:
        for minClusterSize in minClusterSizes:
            size += 1
            print(f'Iteration {size} of {len(minBagSizes)*len(minClusterSizes)}')
            runSRIQ(data=data, minBagSize=minBagSize, minClusterSize = minClusterSize)


def get2Clusters(clusterpaths = None, datapath = '../software/data/splitData.txt', sets = False, title = None, bars = False):

    clusters = [[],[]]
    if clusterpaths is None: clusterpaths = [cp ,cp2]
    for counter, clusterpath in enumerate(clusterpaths):
        clusterfiles = sorted([f for f in listdir(clusterpath) if isfile(join(clusterpath, f)) and f[0] != '.'])
        for file in clusterfiles:
            clusters[counter].append([line.strip() for line in open(f'{clusterpath}/{file}')][1:])

    rows = [inner for outer in clusters[0] for inner in outer]
    columns= [inner for outer in clusters[1] for inner in outer]
    settet = list(set(columns) & set(rows))
    if sets:
        rows = [x for x in rows if x in settet]
        columns = [x for x in columns if x in settet]

    
    newCluster = list()
    for big_cluster in clusters:
        tempList = list()
        for cluster in big_cluster:
            if sets:
                tempList.append([x for x in cluster if x in settet])
            else:
                tempList.append(cluster)
        newCluster.append(tempList)
    clusters1, clusters2 = newCluster[0], newCluster[1]
    
    if bars is True: compareBars(clusters1, clusters2)
    elif bars == 'simi': return (similarity(clusters1, clusters2))
    else: 
        consensusPlot(clusters1, clusters2, title = title)
        return(clusters1, clusters2)







def removeCluster(datapath = None, clusterpath = None):
    df = pd.read_csv(datapath, sep = '\t')
    
    columns = df.columns.tolist()
    clusterfiles = sorted([f for f in listdir(clusterpath) if isfile(join(clusterpath, f)) if f[0] != '.'])
    for counter, file in enumerate(clusterfiles):
        print(f'Run {counter+1}')
        cluster = [line.strip() for line in open(f'{clusterpath}/{file}')][1:]
        fColumns = list(set(columns)- set(cluster))
        fColumns.remove('Gene')
        fColumns.insert(0, 'Gene')
        df[fColumns].to_csv(f'data/expressionData/noCluster{counter+1}.txt', sep = '\t', index = False)
        runSRIQ(data = f'noCluster{counter+1}', studyName = f'No_Cluster_{counter+1}')


# removeCluster()


def compareConsensus(csv = None, cp = None, sets = True,cNum = 5, figSize = (10,10), title = 'test'):
    if csv is None or cp is None:
        print('Set csv and cp first')
    else:
        df = pd.read_csv(csv, sep = ' ')
        clusters1 = list()
        clusters2 = list()
        for i in range(1,cNum +1):
            clusters2.append(df['Assay'][df[f'ConsensusCluster{cNum}'] == i].tolist())
        
        clusterfiles = sorted([f for f in listdir(cp) if isfile(join(cp, f))])
        for file in clusterfiles:
            clusters1.append([line.strip() for line in open(f'{cp}/{file}')][1:])
            
        for counter, cluster in enumerate(clusters2):
            clusters2[counter] = [x.replace('.' , '-' ) for x in cluster]
        rows = [inner for outer in clusters1 for inner in outer]
        columns= [inner for outer in clusters2 for inner in outer]
        settet = list(set(columns) & set(rows))
        if sets:
            clusters1 = [[x for x in cluster if x in settet] for cluster in clusters1]
            clusters2 = [[x for x in cluster if x in settet] for cluster in clusters2]
            
        g = consensusPlot(clusters1, clusters2, title = title, figSize=figSize)
        #return(clusters2, clusters1)

def consensusPlot(clusters1, clusters2, title = 'test', figSize = (7,7)):
    c = np.empty([len(clusters1),len(clusters2)])
    sizes = list()
    for counter, cluster1 in enumerate(clusters1):
        t = list()
        for i, cluster2 in enumerate(clusters2):     
            if len(cluster1) != 0: 
                t.append((len(set(cluster1) & set(cluster2))/len(set(cluster1+cluster2))*100))
                sizes.append(len(set(cluster1+cluster2)))
            else: 
                t.append(0)
                sizes.append(len(set(cluster1+cluster2)))
        c[counter] = t
    xticks, yticks = [x for x in range(1, (c.shape[1]+1))], [x for x in range(1, (c.shape[0]+1))]
    f, ax = plt.subplots(figsize=figSize)
    g = sns.heatmap(c, annot = True, vmax= 100, xticklabels = xticks, yticklabels = yticks)
    g.set(xlabel='SRIQ', ylabel='SRIQ')
    if title is not None: g.set_title(title)
    for counter, t in enumerate(g.texts): 
        if t.get_text() != '0': 
            text = f'{t.get_text()}% [{sizes[counter]}]'
            t.set_text(text)
    
# clusters1, clusters2 = compareConsensus(sets=True)

def similarity(clusters1, clusters2):
    ans = 0
    allClusters = list()
    [[allClusters.append(x) for x in lista] for lista in clusters2]
    [[allClusters.append(x) for x in lista] for lista in clusters1]
    for cluster1 in clusters1:
        for cluster2 in clusters2:
            ans += len(set(cluster1)&set(cluster2))
    return ans/len(set(allClusters))




def r2(cp = (None,)):
    # print (cp[1])
    clusters = list()
    
    if cp[1]:
        #print(cp[0])
        df = pd.read_csv(cp[0], sep = ' ', engine = 'python')
        for i in range(1,5):
            clusters.append(df['Assay'][df['ConsensusCluster4'] == i].tolist())
    else: 
        clusterfiles = sorted([f for f in listdir(cp[0]) if isfile(join(cp[0], f)) if f[0]!='.'])
        for file in clusterfiles:
            clusters.append([line.strip() for line in open(f'{cp[0]}/{file}')][1:])
    collections = list()
    path = f''
    for cluster in clusters:
        cluster = sorted(cluster)
        collections.append(list(combinations(cluster, 2)))
    newCollections = list()
    [[newCollections.append(x) for x in collection] for collection in collections]
    return newCollections

def r1(cp = (None,)):
    # print (cp[1])
    clusters = list()
    
    if cp[1]:
       # print(cp[0])
        df = pd.read_csv(cp[0], sep = ' ', engine = 'python')
        for i in range(1,5):
            clusters.append(df['Assay'][df['ConsensusCluster4'] == i].tolist())
    else: 
        clusterfiles = sorted([f for f in listdir(cp[0]) if isfile(join(cp[0], f)) if f[0]!='.'])
        for file in clusterfiles:
            clusters.append([line.strip() for line in open(f'{cp[0]}/{file}')][1:])
    cluster = list()
    [[cluster.append(x) for x in lista] for lista in clusters]
    return cluster


def r22(mainPath = None, cps = [None], lista = False, reso = None):
    mainC1 = r2(mainPath)
    mainC2 = r1(mainPath) 
    res1 = list()
    res2 = list()
    for cp in cps:
        #print(cp)
        #print((r1(cp)))
        #print('mainC:', len(mainC),'cp:', len(set(r2(cp))), 'mainC&cp:', len(set(mainC)&set(r2(cp))), 'R=',len(set(r2(cp))&set(mainC))/len(mainC))
        res1.append(len(set(r2(cp))&set(mainC1))/len(set(set(mainC1).union(set(r2(cp))))))
        res2.append(len(set(r1(cp))&set(mainC2))/len(set(mainC2).union(set(r1(cp)))))
    return res1, res2


def plotR2(res, res2, names, title = 'test', xTitle = '', yTitle = '', line1 = '', line2 = ''):
    df = pd.DataFrame({'pairwise':res, 'samples':res2}, index = names)
    ax = sns.lineplot(data=df)
    ax.set(ylim= (0,1.1))
    ax.set_title(title)
    ax.set(xticks = names)
    ax.set(xlabel = xTitle,ylabel = yTitle)
    ax.invert_xaxis()
    #plt.show()

def runRobustness(folderPath, title = 'test', xTitle = '', yTitle = '', lista = False, sets = False, returnRes = False):
    folderList = sorted([int(x) for x in os.listdir(folderPath) if x[0] != '.'], reverse = True)
    names = list()
    [names.append(x) for x in folderList]
    names.pop(0)
    folderList = [(folderPath+'/' +str(x)+'/', False) for x in folderList]
    mp = folderList.pop(0)
    #print (mp, folderList)
    res1, res2 = r22(mp, cps = folderList, lista = False)
    if returnRes:
        return res1, res2
    else:
        plotR2(res = res1, res2 = res2, names = names, title = title, xTitle = '', yTitle = '', line1 = 'Common Pairs', line2 = 'Common Samples')

#runRobustness('/Users/jacobkarlstrom/projekt/SRIQ/notebook/data/RobustnessData/BagVariation/3C')

