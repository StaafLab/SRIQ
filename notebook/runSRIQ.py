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
from sklearn.model_selection import train_test_split
import pandas as pd
import numpy as np
from os import listdir
from os.path import isfile, join
import seaborn as sns
np.random.seed(10)

def runSRIQ(data, cutOff = None, permutations = 10000, iterations = 10, minBagSize=500, minClusterSize = 0):
    #data: string with the directory to the csv file which you want to run
    #cutOff: list of the cutoffs
    output = ''
    resources = '../software/VRLA/resources/test.properties'
    if cutOff is None: cutOff = [0.9,0.89,0.88,0.87,0.86,0.85,0.84,0.83,0.82,0.81,0.80,0.79,0.78,0.77,0.76,0.75,0.74,0.73,0.72,0.71,0.7,0.69,0.68,0.67,0.66,0.65,0.64,0.63,0.62,0.61,0.6,0.59,0.58,0.57,0.56,0.55,0.54,0.53,0.52,0.51,0.5,0.49,0.48,0.47,0.46,0.45,0.44,0.43,0.42,0.41,0.4,0.39,0.38,0.37,0.36,0.35,0.34]
    with open(resources) as file:
        for line in file.readlines():
            if 'inFileName' in line:
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
    bashCommand = "cd ../software/VRLA && java -jar -Xmx4g VRLA.jar"
    os.system(bashCommand)
    print('Done!')
    
# runSRIQ(data = 'test_mc_log2var(80)', minBagSize = 1500)
def splitRunSRIQ(data, frac, per):
    data = f'../software/data/{data}'
    df = pd.read_csv(data, sep = '\t')
    remove_n = int(df.shape[1]*(frac))
    col = df['Gene']
    columns = np.random.choice(df.columns, remove_n, replace=False)
    columns[0] = 'Gene'
    df = df[columns]
    print(df.shape)
    df.to_csv(f'../software/data/splitData.txt', sep = '\t', index = False)
    runSRIQ('splitData', permutations = per)

# lista = [(0.9,9998),(0.8,9997),(0.7,9996),(0.6,9995)]
# for ele in lista:
#     splitRunSRIQ(data = 'test_mc_log2var(80).txt', frac = ele[0], per=ele[1])

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
    # cp = '/home/jacob/examensarbete/software/output/allt/sparad/10000/QC_Spiral(false)/dist(0.43)/4'
    # cp2 = '/home/jacob/examensarbete/software/output/allt/sparad1/10000/QC_Spiral(false)/dist(0.45)/4'
    cp = '/home/jacob/examensarbete/software/output/VRLA_test_10000itr_1400var_10r/10000/QC_Spiral(false)/dist(0.47)/4'
    # cp2 = '/home/jacob/examensarbete/software/output/VRLA_test_10000itr_1500var_10r/10000/QC_Spiral(false)/dist(0.47)/4'
    cp2 = '/home/jacob/examensarbete/software/output/VRLA_test_10000itr_900var_10r/10000/QC_Spiral(false)/dist(0.45)/4'
    if clusterpaths is None: clusterpaths = [cp ,cp2]
    for counter, clusterpath in enumerate(clusterpaths):
        clusterfiles = sorted([f for f in listdir(clusterpath) if isfile(join(clusterpath, f))])
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
    
    if bars: compareBars(clusters1, clusters2)
    else: consensusPlot(clusters1, clusters2, title = 'SRIQ 5 clusters vs 4 cluster')
    return clusters1, clusters2

cp = '/home/jacob/examensarbete/software/output/VRLA_test_10000itr_1200var_10r/10000/QC_Spiral(false)/dist(0.48)/5'
# cp2= '/home/jacob/examensarbete/software/output/allt/sparad/10000/QC_Spiral(false)/dist(0.43)/4'
cp2 = '/home/jacob/examensarbete/software/output/VRLA_test_10000itr_1200var_10r/10000/QC_Spiral(false)/dist(0.45)/4'
get2Clusters([cp, cp2],sets= False, bars = False)



def consensusPlot(clusters1, clusters2, title = 'test'):
    c = np.empty([len(clusters1),len(clusters2)])
    
    for counter, cluster1 in enumerate(clusters1):
        t = list()
        for i, cluster2 in enumerate(clusters2):     
            if len(cluster1) != 0: t.append(len(set(cluster1) & set(cluster2))/len(set(cluster1+cluster2))*100)
            else: t.append(0)
        c[counter] = t
    xticks, yticks = [x for x in range(1, (c.shape[1]+1))], [x for x in range(1, (c.shape[0]+1))]
    g = sns.heatmap(c, annot = True, vmax= 100, xticklabels = xticks, yticklabels = yticks)
    g.set(xlabel='4 clusters', ylabel='5 clusters')
    if title is not None: g.set_title(title)
    for t in g.texts: 
        if t.get_text() != '0': 
            t.set_text(t.get_text() + " %")


def removeCluster(datapath = None, clusterpath = None):
    if clusterpath == None: clusterpath = '/home/jacob/examensarbete/software/output/allt/sparad/10000/QC_Spiral(false)/dist(0.45)/4'
    if datapath == None: datapath = '/home/jacob/examensarbete/software/data/test_mc_log2var(80).txt'
    df = pd.read_csv(datapath, sep = '\t')
    
    columns = df.columns.tolist()
    clusterfiles = sorted([f for f in listdir(clusterpath) if isfile(join(clusterpath, f))])
    for counter, file in enumerate(clusterfiles):
        print(f'Run {counter+1}')
        cluster = [line.strip() for line in open(f'{clusterpath}/{file}')][1:]
        fColumns = list(set(columns)- set(cluster))
        fColumns.remove('Gene')
        fColumns.insert(0, 'Gene')
        df[fColumns].to_csv(f'../software/data/noCluster{counter+1}.txt', sep = '\t', index = False)
        runSRIQ(data = f'noCluster{counter+1}', permutations = 10000+counter)


# removeCluster()


def compareConsensus(csv = 'data/Test_Pearson_WardD_ConsensusClasses2.txt', cp = '/home/jacob/examensarbete/software/output/VRLA_test_10000itr_1200var_10r/10000/QC_Spiral(false)/dist(0.71)/4', sets = True):
    df = pd.read_csv(csv, sep = ' ')
    clusters1 = list()
    clusters2 = list()
    for i in range(1,5):
        clusters2.append(df['Assay'][df['ConsensusCluster4'] == i].tolist())
    
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
        
    g = consensusPlot(clusters1, clusters2, title = 'Sets of SRIQ vs consensus NC = 4')
    return(clusters2, clusters1)
    
# clusters1, clusters2 = compareConsensus(sets=True)

def compareBars(clusters1, clusters2):
    c = [len(x) for x in clusters1]
    [c.append(len(x)) for x in clusters2]
    c2 = ['SRIQ%' if x/(len(c)-1) < 0.5 else 'Consensus' for x in range(len(c))]
    c3 = [x for x in range(1,5)]
    [c3.append(x) for x in range(1,5)]
    df = pd.DataFrame({'Cluster sizes':c,'run':c2, 'Cluster':c3})
    
    ax = sns.barplot(data = df, x='Cluster', y='Cluster sizes', hue='run')

# compareBars(clusters1, clusters2)



def r2(cp = (None,)):
    # print (cp[1])
    clusters = list()
    
    if cp[1]:
        print(cp[0])
        df = pd.read_csv(cp[0], sep = ' ', engine = 'python')
        for i in range(1,5):
            clusters.append(df['Assay'][df['ConsensusCluster4'] == i].tolist())
    else: 
        clusterfiles = sorted([f for f in listdir(cp[0]) if isfile(join(cp[0], f))])
        for file in clusterfiles:
            clusters.append([line.strip() for line in open(f'{cp[0]}/{file}')][1:])
    collections = list()
    path = f''
    for cluster in clusters:
        collections.append(list(combinations(cluster, 2)))
    newCollections = list()
    [[newCollections.append(x) for x in collection] for collection in collections]
    return newCollections


def r22(mainPath = None, cps = [None]):
    mainC = r2(mainPath)
    
    res = list()
    for cp in cps:
        # print(cp)
        print('mainC:', len(mainC),'cp:', len(set(r2(cp))), 'mainC&cp:', len(set(mainC)&set(r2(cp))), 'R=',len(set(r2(cp))&set(mainC))/len(mainC))
        res.append(len(set(r2(cp))&set(mainC))/len(set(mainC).union(set(r2(cp)))))
    return res


def plotR2(res):
    df = pd.DataFrame({'Similarity':res, 'Bag size':['9000', '12000', '18000']})
    ax = sns.lineplot(data=df, x="Bag size", y="Similarity")
    ax.set(ylim= (0,1.1))
    ax.set_title('Pair similarity 5 clusters: 1500 vs:')

# mp = ('/home/jacob/examensarbete/software/output/VRLA_test_10000itr_1500var_10r/10000/QC_Spiral(false)/dist(0.49)/5', False)
# cp2 = ('data/consensusData/Test_Pearson_WardD_ConsensusClasses.txt', True)
# cp3 = ('data/consensusData/Test_Pearson_WardD_ConsensusClasses2.txt', True)
# cp4 = ('data/consensusData/Test_Pearson_WardD_ConsensusClasses3.txt', True)
# # cp2 = ('/home/jacob/examensarbete/software/output/VRLA_test_10000itr_1400var_10r/10000/QC_Spiral(false)/dist(0.5)/5', False)
# mp = ('/home/jacob/examensarbete/software/output/VRLA_test_10000itr_1200var_10r/10000/QC_Spiral(false)/dist(0.48)/5', False)
# # cp4 = ('/home/jacob/examensarbete/software/output/VRLA_test_10000itr_900var_10r/10000/QC_Spiral(false)/dist(0.47)/5', False)
# # cp5 = ('/home/jacob/examensarbete/software/output/VRLA_test_10000itr_406var_10r/10000/QC_Spiral(false)/dist(0.47)/5', False)
# # cp6 = ('/home/jacob/examensarbete/software/output/VRLA_test_10000itr_400var_10r/10000/QC_Spiral(false)/dist(0.44)/5', False)
# cps =[cp2, cp3, cp4]

# # cps = ['/home/jacob/examensarbete/software/output/VRLA_test_10000itr_1500var_10r/10000/QC_Spiral(false)/dist(0.47)/4', '/home/jacob/examensarbete/software/output/VRLA_test_10000itr_1400var_10r/10000/QC_Spiral(false)/dist(0.47)/4', '/home/jacob/examensarbete/software/output/VRLA_test_10000itr_1200var_10r/10000/QC_Spiral(false)/dist(0.45)/4', '/home/jacob/examensarbete/software/output/VRLA_test_10000itr_900var_10r/10000/QC_Spiral(false)/dist(0.45)/4','/home/jacob/examensarbete/software/output/VRLA_test_10000itr_406var_10r/10000/QC_Spiral(false)/dist(0.43)/4', '//home/jacob/examensarbete/software/output/allt/test/10000/QC_Spiral(false)/dist(0.43)/4']
# # cps = [(cp, False) for cp in cps]
# # mp = cps[0]
# res = r22(mp, cps = cps)
# plotR2(res)
