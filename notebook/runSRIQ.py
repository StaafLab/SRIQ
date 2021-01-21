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
    

def splitRunSRIQ(data, jupyter = True):
    data = f'../software/data/{data}'
    df = pd.read_csv(data, sep = '\t')
    remove_n = int(df.shape[1] / 2)
    col = df['Gene']
    columns = np.random.choice(df.columns, remove_n, replace=False)
    columns[0] = 'Gene'
    df = df[columns]
    df.to_csv(f'../software/data/splitData.txt', sep = '\t', index = False)
    runSRIQ('splitData')



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


def consensusPlot(clusterpaths = None, datapath = '../software/data/splitData.txt'):

    clusters = [[],[]]
    if clusterpaths == None: clusterpaths = ['/home/jacob/examensarbete/software/output/allt/sparad1/10000/QC_Spiral(false)/dist(0.45)/4','/home/jacob/examensarbete/software/output/allt/sparad/10000/QC_Spiral(false)/dist(0.45)/4']
    for counter, clusterpath in enumerate(clusterpaths):
        clusterfiles = sorted([f for f in listdir(clusterpath) if isfile(join(clusterpath, f))])
        for file in clusterfiles:
            clusters[counter].append([line.strip() for line in open(f'{clusterpath}/{file}')][1:])

    rows = [inner for outer in clusters[0] for inner in outer]
    columns= [inner for outer in clusters[1] for inner in outer]
    settet = list(set(columns) & set(rows))
    rows = [x for x in rows if x in settet]
    columns = [x for x in columns if x in settet]
    data = list()
    for row in rows:
        res = list()
        flag = False
        for column in columns:
            for i in range(len(clusters[0])):
                if row in clusters[0][i] and row in clusters[1][i] and column in clusters[0][i] and column in clusters[1][i]:
                    flag = True
            if flag == True: res.append(1)
            else: res.append(0)
            flag = False
        data.append(res)
    df = pd.DataFrame(data, columns = columns, index = rows)
    return df

sns.heatmap(consensusPlot())

def removeCluster(datapath = None, clusterpath = None):
    if clusterpath == None: clusterpath = '/home/jacob/examensarbete/software/output/sparad/10000/QC_Spiral(false)/dist(0.45)/4'
    if datapath == None: datapath = '/home/jacob/examensarbete/software/data/test_mc_log2var(80).txt'
    df = pd.read_csv(datapath, sep = '\t')
    
    columns = df.columns.tolist()
    clusterfiles = sorted([f for f in listdir(clusterpath) if isfile(join(clusterpath, f))])
    for counter, file in enumerate(clusterfiles):
        cluster = [line.strip() for line in open(f'{clusterpath}/{file}')][1:]
        fColumns = list(set(columns)- set(cluster))
        fColumns.remove('Gene')
        fColumns.insert(0, 'Gene')
        df[fColumns].to_csv(f'../software/data/noCluster{counter+1}.txt', sep = '\t', index = False)
        runSRIQ(data = f'test_mc_log2', permutations = 10000+counter)

        