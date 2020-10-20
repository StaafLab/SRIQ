#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 10:04:27 2020

@author: jacob
"""
import numpy as np

dataDirectory = '/home/jacob/examensarbete/software/data/LUAD_GeneExpData_mc_log2var(50).txt'
clusterDirectory = '/home/jacob/examensarbete/R/Data/5/5_5C_Samples.txt'

geneList = []
finalGeneList = list()


with open(clusterDirectory) as file:
    for line in file.readlines():
        line = line.replace('\n' ,'')
        geneList.append(line)
    
    
    
with open(dataDirectory) as file:
    genes = file.readline()
    genes = genes.split()
    
    for counter, gene in enumerate(genes):
        if gene in geneList:
            finalGeneList.append(counter)
    
    for line in file.readlines():
        line = line.split()
        for counter, ele in enumerate(line):
            if counter in finalGeneList or counter == 0:
                if counter == 0: geneList.append(ele)
                else: geneList.append(float(ele))

arr = np.array(geneList)
arr = arr.reshape(10046,len(finalGeneList))

np.savetxt("cluster5Data.csv", arr, delimiter=",",fmt ='%s')