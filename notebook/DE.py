#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 09:04:10 2020

@author: jacob
"""
import mygene
import json
import requests
import os
import sys
import math
from numpy import *
import numpy as np
import pandas as pd
import seaborn as sns
from umap import UMAP
import matplotlib.pyplot as plt
from bioinfokit.analys import get_data, stat
from sklearn.decomposition import PCA
from scipy.stats import ttest_ind
import statsmodels as sm
from tqdm.auto import tqdm
from scipy.stats import normaltest
from scipy.stats import mannwhitneyu
from lifelines import CoxPHFitter
from lifelines import KaplanMeierFitter as kmf
from pandas.api.types import is_numeric_dtype



class networkAnalysis():
    def __init__(self, newDf = None):
        #Initiates the working dataframe for manual initiation, not recomended.
        if newDf:
            self.df = newDf
            self.tDf = self.df.transpose()
        self.eList = None
        self.pDf = None
        self.lfcDf = None
        self.rDf = None
        self.symbolDf = None
        self.clinicalDf = None
        self.metaDf = None
        self.eDf = None
        
        
    def renameCol(self, oldCol, newCol):
        self.df = self.df.rename(columns={oldCol:newCol})
        self.tDf = self.df.transpose()
        
    def log2med(self, i):
        self.log2(i, 'l2m')
        self.medianCenter(i, 'l2m')
    
    def medlog2(self, i):
        self.medianCenter(i, 'ml2')
        self.log2(i,',ml2')
        
    def toCsv(self, path ,sep = ','):
        self.df.to_csv(path, sep = sep, index = False)
  
    def readCsv(self, csvpath = None, readType = 'csv', newIndex = False, sep = ',', fpkm = False):
            if readType == 'csv':
                self.df = pd.read_csv(csvpath, sep = sep)
                self.tDf = self.df.transpose()  
            elif readType == 'fpkm':
                self.fpkmDf = pd.read_csv(csvpath, sep = sep).set_index('Unnamed: 0').transpose()
                meanList = list()
                for gene in list(self.fpkmDf):
                    meanList.append(self.fpkmDf[gene].mean())
                self.fpkmDf = self.fpkmDf.transpose()
                self.fpkmDf['FPKM'] = meanList
                self.fpkmDf = self.fpkmDf['FPKM']
            elif readType == 'centroid':
                self.cDf = pd.read_csv(csvpath, sep = sep, index_col = 'Unnamed: 0')
            elif readType=='lfc':
                self.rawDf = pd.read_csv(csvpath, sep = sep, index_col = 'Unnamed: 0')
            elif readType == 'meta-gene':
                self.metaDf = pd.read_excel(r'data/clinicalData/130932_1_supp_2656123_nbr7l8.xlsx')
                self.metaDf = self.metaDf.iloc[:, 0:2]
            elif readType == 'Clinical':
                if csvpath is None:
                    self.clinicalDf = pd.read_csv(r'data/clinicalData/clinical_PANCAN_patient_with_followup.tsv', sep = '\t', engine = 'python', index_col = 'bcr_patient_barcode')
                else:
                    self.clinicalDf = pd.read_csv(csvpath, sep = '\t', engine = 'python', index_col = 'bcr_patient_barcode')
            if newIndex:
                self.df = self.df.set_index(newIndex)
                self.tDf = self.df.transpose()
                
            
                

    def test(self, i):
        #preprocessing of data with filtering of fpkm < 1 built in
        self.df.iloc[:,1:] = self.df.iloc[:,1:].where(self.df.iloc[:,1:] >= 1, 1)
        self.tDf = self.df.transpose()
        self.tDf.iloc[i:] = self.tDf.iloc[i:].transform(lambda x: x/x.median() if x.median() != 0 else x)
        self.tDf.iloc[i:] = np.log2(self.tDf.iloc[i:].astype('float64'))
        self.df = self.tDf.transpose()
        
        
        
    def medianCenter(self, i, method = 'l2m'):
        if method == 'l2m': self.tDf.iloc[i:] = self.tDf.iloc[i:].transform(lambda x: x-x.median())
        else: self.tDf.iloc[i:] = self.tDf.iloc[i:].transform(lambda x: x/x.median() if x.median() != 0 else x)
        self.df = self.tDf.transpose()
        
    def log2(self, i, method = 'l2m'):
        if method == 'l2m': self.df.iloc[:,i:] = np.log2(self.df.iloc[:,i:].transform(lambda x: x+1))
        else : self.df.iloc[:,i:] =  np.log2(self.df.iloc[:,i:].replace(0, np.nan)).replace(np.nan, 0)
        self.tDf = self.df.transpose()
        
        
        
    def filterVariantGenes(self, start = 0, top = 1, bottom = 0):
        #Filter the least variant genes based on threshold ranging 0 to 1
        self.filterDf = self.df.iloc[:,start:].loc[(self.df.iloc[:,start:].var(axis = 1) < self.df.iloc[:,start:].var(axis = 1).quantile(top)) & (self.df.iloc[:,start:].var(axis = 1) > self.df.iloc[:,start:].var(axis = 1).quantile(bottom))]
    
    def readSRIQ(self, csvpath, clusterpath, columnname = 'gene_id'):
        
        self.df = pd.read_csv(csvpath, delimiter = '\t', index_col = columnname)
        self.sortedClusters = list()
        self.allClusters = list()
        for file in sorted(os.listdir(clusterpath)):
            lista = [line.strip() for line in open(f'{clusterpath}/{file}')][1:]
            self.allClusters += lista
            self.sortedClusters.append(lista)
        self.df = self.df[self.allClusters]
        clusterRow = [inner for outer in [[counter+1]*len(x) for counter, x in enumerate(self.sortedClusters)] for inner in outer]
        clusterRow = pd.Series(data = clusterRow, name = 'Clusters', index = self.df.columns)
        self.df = self.df.append(clusterRow)
        self.df = self.df[(self.df.T != 0).any()]
        self.tDf = self.df.transpose()
        lut = dict(zip(self.tDf['Clusters'].unique(), 'bgrcmykw'))
        self.col_colors = self.tDf['Clusters'].map(lut)
        ax = sns.countplot(self.tDf['Clusters'])
        self.filterDf = self.df
        self.clusterNum = clusterpath[-1]
        f = lambda l: ['-'.join(s.split('-')[0:3]) for s in l]
        self.samples = f(self.df.columns.tolist())
            
               
    def setlog2Df(self, df):
        self.df = df
        self.tDf = self.df.transpose()
        

    def screePlot(self, Clustered = True, demo = False):
        #Calculates the pca model and plots the scree plot
        pca = PCA()
        if Clustered: 
            pca.fit(self.tDf.drop(['Clusters'],axis=1))
            self.pca_data = pca.transform(self.tDf.drop(['Clusters'],axis=1))
        else: 
            pca.fit(self.tDf)
            self.pca_data = pca.transform(self.tDf)
        
        per_var = np.round(pca.explained_variance_ratio_* 100, decimals=1)
        self.labels = ['PC' + str(x) for x in range(1, len(per_var)+1)]
        #If its the demorun it will save the pca model into a dataframe for later use
        if demo:
            numList = [i+1 for i in range(len(per_var))]
            lista = [demo]*len(per_var)
            self.testDf = pd.DataFrame({'per_var':per_var, 'cutoff':lista, 'num':numList})
            self.sDf = self.sDf.append(self.testDf, ignore_index=True)
        #If it isnt the demorun it will plot the scree plot
        else:
            plt.bar(x=range(1,len(per_var)+1), height=per_var)
            plt.ylabel('Percentage of Explained Variance')
            plt.xlabel('Principal Component')
            plt.title('Scree Plot')
            plt.show()
        


        
    def pcaPlot(self, Clustered = True):
        #Scatterplots the pca model
        self.pca_df = pd.DataFrame(self.pca_data,  columns=self.labels)
        if Clustered:
            self.pca_df = self.pca_df.join(self.tDf.filter(['Clusters']).reset_index())
            sns.relplot(x="PC1", y="PC2", data=self.pca_df, hue='Clusters', style='Clusters');
        else:
            sns.relplot(x="PC1", y="PC2", data=self.pca_df);
        plt.show()
        
        
    def Umap(self,n, demo = False):
        #Does a umap model based on the pca model, n = number of desired principal components
        X = self.tDf.values[:,0:(self.tDf.shape[1]-1)]
        Y = self.tDf.values[:,self.tDf.shape[1]-1]
        X_reduced = PCA(n_components = 30).fit_transform(X)
        print("Performing Uniform Manifold Approximation and Projection (UMAP) ...")
        model = UMAP(n_neighbors = 30, min_dist = 0.3, n_components =n )
        self.umap = model.fit_transform(X_reduced)
        if demo:
            cutoff = [demo]*len(self.umap[:,0])
            self.tempDf = pd.DataFrame({'comp 1':self.umap[:,0],'comp 2':self.umap[:,1], 'cutoff':cutoff ,'hue':Y})
            self.uDf = self.uDf.append(self.tempDf,ignore_index = True)
        else:
            sns.scatterplot(x=self.umap[:,0], y=self.umap[:,1], palette = 'viridis', size = 1, style = Y, hue = Y)

        
    def transformDf(self, cpm = False):
        #Old method for un log transforming to calculate cpm for filtering.
        self.transDf = self.df.transform(lambda x: 2**x)
        self.tTransDf = self.transDf.transpose()
        if cpm:
            lista = list()
            for sample in self.allClusters:
                lista.append(self.transDf[sample].sum()/1000000)
            self.tTransDf['CPM'] = lista
            self.newDf = self.tTransDf.iloc[:,0:-1].div(self.tTransDf['CPM'], axis = 0)

    
    def diffGeneAnalysis(self, transformed = True, test = 'non-parametric'):
        
        #Tests each gene for its significance in the different clusters
        self.pDf = pd.DataFrame(index = [x.split('.')[0] for x in self.filterDf.index.tolist()])
        if test == 't-test':
            for counter, columns in enumerate(tqdm(self.sortedClusters)):
                res = list()
                set2 = set(self.filterDf) - set(columns)
                df1, df2 = self.filterDf[columns], self.filterDf[set2]
                for i in tqdm(range(df1.shape[0])):
                    res.append(ttest_ind(df1.iloc[i].values, df2.iloc[i].values,equal_var = False)[1])
                    if np.isnan(res[-1]):
                        res.pop(-1)
                        res.append(1)
                self.pDf[f'Cluster {counter +1}'] = res
        elif test == 'mannwhitneyu':
            for counter, columns in enumerate(tqdm(self.sortedClusters)):
                res = list()
                set2 = set(self.filterDf) - set(columns)
                df1, df2 = self.filterDf[columns], self.filterDf[set2]
                for i in tqdm(range(df1.shape[0])):
                    res.append(mannwhitneyu(df1.iloc[i].values, df2.iloc[i].values)[1])
                self.pDf[f'Cluster {counter +1}'] = res
        elif test == 'kruskal':
            for counter, columns in enumerate(tqdm(self.sortedClusters)):
                res = list()
                set2 = set(self.filterDf) - set(columns)
                df1, df2 = self.filterDf[columns], self.filterDf[set2]
                for i in tqdm(range(df1.shape[0])):
                    res.append(mannwhitneyu(df1.iloc[i].values, df2.iloc[i].values)[1])
                self.pDf[f'Cluster {counter +1}'] = res
        #CORRECTED
        print('Performing FDR correction')
        for i in range(len(self.sortedClusters)):
            self.pDf['Cluster {} corr'.format(i+1)] = sm.stats.multitest.fdrcorrection(self.pDf['Cluster {}'.format(i+1)].array, alpha=0.05, method='indep', is_sorted=False)[1]
        print('Filtering out genes')
        self.pDf = self.pDf.loc[(self.pDf.iloc[:,4:] < 0.01).sum(axis = 1) < len(self.sortedClusters)]
        bla = pd.DataFrame(columns = [f'Cluster {i+1} corr:' for i in range(len(self.sortedClusters))])
        # Sets the lesser significant genes to a p value of 1
        tempDf = self.pDf.iloc[:,len(self.sortedClusters):].min(axis = 1)
        for i in tqdm(range(self.pDf.shape[0])):
            bla = bla.append(self.pDf.iloc[:,len(self.sortedClusters):].iloc[i].transform(lambda x: x if x == tempDf.iloc[i] else x/x))
        self.pDf.iloc[:,len(self.sortedClusters):] = bla
        
            
    def filtering(self, threshold = 2, filteringType = True, start = 1 ,transformed = True, csvpath = 'data/test.csv', sep = '.'):
        #Filters based on fpkm value. Not needed if using the test metho use threshold = 1
        if self.pDf is not None:
            if filteringType == 'fpkm': self.readCsv(csvpath, readType = 'fpkm',  sep = ',', fpkm = True)
            elif filteringType == 'log2fold': self.logFoldChange()
            #
            genes = list()
            #Fethes the significant genes
            for i in range(len(self.sortedClusters)):
                genes.append(self.pDf[self.pDf['Cluster {} corr'.format(i+1)] < 0.01].index.tolist())
            self.tList = list()
            #Filters the genes based of fpkm
            if filteringType == 'fpkm':
                for i in range(len(self.sortedClusters)):
                    finalList = genes[i]
                    lista = list()
                    for gene in finalList:
                        if  (self.fpkmDf[gene] > threshold):
                            lista.append(gene)
                    self.tList.append(lista)
            if filteringType == 'log2fold':
                for i in range(len(self.sortedClusters)):
                    finalList = genes[i]
                    #I dont fucking know
                    self.tList.append(list(set(genes[i]) & set(list(self.lfcDf.loc[(self.lfcDf[f'Cluster {i+1}']) > threshold].index))))
            for counter, genes in enumerate(self.tList):
                print(f'Cluster {counter+1} {len(genes)}')
                
            #Produces an enrichment list which can be used for enrichR
            self.eList = list()
            sep = sep
            for i, genes in enumerate(self.tList):
                lista = list()
                for j, gene in enumerate(genes):
                    lista.append(gene.split(sep,1)[0])
                self.eList.append(lista)
        else: print('Error: Run clusterTTest first')
            
    def log2filtering(self):
        #Filters groups with log2fold change lesser than 2
        self.df = self.df.filter(items = list(self.lfcDf[abs(self.lfcDf > 2)].dropna(axis=0,how='all').index), axis = 0)
        self.tDf = self.df.transpose()

            
    def clusterMap(self, cluster = 'all', clustered = True, vmin = -1, vmax = 1, row_cluster = False):
        #plots a clustermap of the working df
        if clustered: 
            if cluster == 'all':
                #Creates the dataframe to be plotted by filtering on the enriched gene set and stacking them
                #on top of each other cluster wise
                self.filterDf.index = [x.split('.')[0] for x in self.filterDf.index.tolist()]
                self.mapDf = pd.DataFrame(columns = list(self.df))
                cols, newList = 'bgrcmykw', list() 
                for i, lista in enumerate(self.tList):
                    self.mapDf = pd.concat([self.mapDf, self.filterDf.filter(items = lista, axis = 'index')])
                    newList.append([cols[i]]*len(lista))
                newList = [inner for outer in newList for inner in outer]
                lista = [inner for outer in self.tList for inner in outer] #turns list of lists into a list
                self.row_colors = pd.DataFrame(index = lista, data = newList)
                pd.DataFrame(index = lista, data = newList)
                self.row_colors = self.row_colors.rename(columns={0:'Clusters'})
                self.row_colors = self.row_colors['Clusters']
                # for i, genes in enumerate(tqdm(self.tList)):
                #     self.mapDf = pd.concat([self.mapDf, self.tDf.filter(items=genes)], axis = 1)
                ax = sns.clustermap(data = self.mapDf, cmap = 'vlag', metric = 'correlation', method = "average",vmax= vmax, vmin = vmin,col_cluster=False, col_colors=self.col_colors, row_colors = self.row_colors, row_cluster = row_cluster)
            else:
                sns.clustermap(data = self.filterDf.filter(items=self.tList[(cluster)], axis = 'index'), cmap = 'vlag', metric = 'correlation', method = "average",vmax= vmax, vmin = vmin,col_cluster=False, col_colors=self.col_colors, row_cluster = row_cluster)
        else: sns.clustermap(data = self.filterDf, cmap = 'vlag', metric = 'correlation', method = "average",vmax= vmax, vmin = vmin,col_cluster=False, col_colors=self.col_colors, row_cluster = row_cluster)
        
    
    def enrichR(self, secondRun = False, dbs = ['Enrichr_Libraries_Most_Popular_Genes']):
        #Fetches enrichment info from enrichR against different databases.
        #Fetches user list id for api calls
        if self.eList is not None:
            c = [0,1,2,3,4,5,6,7,8,'db', 'cluster']
            self.rDf = pd.DataFrame(columns = c)
            for counter, cluster in enumerate(self.eList):
                print(f'Running for cluster {counter+1}')
                ENRICHR_URL = 'http://maayanlab.cloud/Enrichr/addList'
                genes_str = '\n'.join(cluster)
                description = 'Example gene list'
                payload = {
                    'list': (None, genes_str),
                    'description': (None, description)
                }
                response = requests.post(ENRICHR_URL, files=payload)
                if not response.ok:
                    raise Exception('Error analyzing gene list')
                user_list_id = json.loads(response.text)['userListId']
                #Fetches values from each database and saved them longformat into rDf
                if secondRun is True: runDbs = dbs[counter]
                else: runDbs = dbs
                for db in runDbs:
                    print(db)
                    ENRICHR_URL = 'http://maayanlab.cloud/Enrichr/enrich'
                    query_string = '?userListId=%s&backgroundType=%s'
                    response = requests.get(
                        ENRICHR_URL + query_string % (user_list_id, db)
                        )
                    if not response.ok:
                        raise Exception('Error fetching enrichment results')
                        
                    data = json.loads(response.text)
                    testDf = pd.DataFrame(data[db])
                        
                    testDf[2] = testDf[2].transform(lambda x: 1/x)
                    testDf[2] = np.log10(testDf[2])
                    testDf = testDf.head()
                    testDf['db'] = [db]*len(testDf.index.tolist())
                    testDf['cluster'] = [str(counter+1)]*len(testDf.index.tolist())
                    self.rDf = self.rDf.append(testDf)
            self.rDf = self.rDf.sort_values(by = [2], ascending = False)
        else: print('Error: Run the filtering module')
            
    def plotBar(self):
        #plots the enirchment result from rDf
        if self.rDf is not None:
            # g = sns.catplot(x= 2, y=1,kind='bar', data = self.rDf, row = 'cluster')
            size = len(self.rDf[self.rDf['cluster'] == '1'].index.tolist())
            fig, axs = plt.subplots(len(self.sortedClusters), figsize = (size,size))
            for i in range(len(self.sortedClusters)):
                sns.barplot(data = self.rDf[self.rDf['cluster'] == str(i+1)], x = 2, y = 1, ax = axs[i]).set_title(f'Cluster {i+1}')
            fig.tight_layout()
        else: print ('Run the enrichR module')
        
    def demoRun(self, cutoffs, csvpath, clusterpath, columnname = 'Gene'):
        #Allows the user to get an overview over the SRIQ output
        #cutoffs: dictionary of cutoff:[cluster values]
        self.sDf = pd.DataFrame({'per_var', 'cutoff','num' })
        self.uDf = pd.DataFrame(columns = ['comp 1','comp 2','cutoff', 'hue'])
        lista = list()
        
        #Fetches the information from the specified clusters
        for i, cutoff in enumerate(cutoffs):
            print(f'Run {i+1} out of {len(cutoffs)}')
            for j, clust in enumerate(cutoffs[cutoff]):
                cpath = clusterpath.split('/')
                cpath[-2] = 'dist({})'.format(cutoff)
                cpath[-1] = '{}'.format(clust)
                cpath = '/'.join(cpath)
                lista.append(cpath[-8:])
                self.readSRIQ(csvpath, cpath, columnname)
                self.screePlot(demo = cpath[-8:])
                self.Umap(3,demo = cpath[-8:])
                
        #plots the information from the fetched info
        g = sns.catplot(x= 'num', y='per_var',kind='bar', data = self.sDf, col='cutoff',col_wrap= 2, sharey = False, sharex = False)
        g2 = sns.relplot(x= 'comp 1', y='comp 2', data= self.uDf, style = 'hue', hue='hue', col='cutoff',col_wrap=2, sharey = False, sharex = False)
    
        
    def ensemble2gene(self):
        #Converts the enrichment list from ensembl to gene symbol
        for i, ele in enumerate(self.eList):    
            self.eList[i] = self.ens2symHelper(ele)['symbol'].dropna().tolist()

    def metaGenes(self):
        if self.symbolDf is None: self.ensembl2symbol()
        if self.metaDf is None: self.readCsv(readType = 'meta-gene')
        df = self.metaDf
        df = df.iloc[:, 0:2]
        self.boxDf = pd.DataFrame(columns = ['Values', 'Cluster', 'Metagene'])
        for x in df['Network cluster'].unique():
            print(x)
            for i, cluster in enumerate(self.sortedClusters):
                values = (self.symbolDf[cluster].filter(items = list(df[df['Network cluster'] == x]['Gene Symbol']), axis = 0).mean(axis = 1).values)
                tempDf = pd.DataFrame({'Values': values, 'Cluster':[i+1 for c in range(len(values))], 'Metagene':[x for c in range(len(values))]})
                self.boxDf = pd.concat([self.boxDf, tempDf])
        ax = sns.catplot(kind = 'box', data = self.boxDf, x = 'Cluster', y = 'Values', col = 'Metagene', col_wrap = 2)
    
        
        
    def ens2symHelper(self, genes):
        headers = {'content-type': 'application/x-www-form-urlencoded'}
        params = f'q={",".join(genes)}&scopes=ensembl.gene&fields=symbol'
        res = requests.post('http://mygene.info/v3/query', data=params, headers=headers)
        data = json.loads(res.text)
        return pd.json_normalize(data).set_index('query')
        
        
    def ensembl2symbol(self):
        print('This may produce duplicate indexes which will be removed, keeping the first occurence')
        newIndex = [x.split('.')[0] for x in list(self.df.index)]
        self.df.index = newIndex
        self.tDf = self.df.transpose()
        self.genes = self.ens2symHelper(newIndex)
        newIndex = (self.genes['symbol'].dropna())
        self.symbolDf = self.df.reindex(list(newIndex.index))
        self.symbolDf.index = list(newIndex)
        self.symbolDf = self.symbolDf[~self.symbolDf.index.duplicated(keep='first')]
        self.tSymbolDf = self.symbolDf.transpose()
        self.tSymbolDf['Clusters'] = self.tDf['Clusters']

        
                
        
    def centroids(self, centroidpath = 'data/wilkerson.2012.LAD.predictor.centroids.csv', ensembl = True, method = 'pearson'):
        #Purpose: calculating centroids on cancer data
        #read the centroid predictor data set
        self.readCsv(csvpath = centroidpath, readType = 'centroid')
        if self.symbolDf is None:
            #If the working data is in ensemble, it will find the matchning gene-ensembldids.
            #Finallist consist of the matching genes
            self.ensembl2symbol()
        finalList = set(self.symbolDf.index) & set(self.cDf.index)
        #ctDf consists of the working dataframe with only centroid predicting genes.
        self.ctDf = self.symbolDf
        self.ctDf = self.ctDf.filter(items = finalList, axis = 'index')
        #Transform cDf to a series of each cancers type's sum
        self.cDf = self.cDf.filter(items = finalList, axis = 'index')
        
        cList = list()
        #Find what centroid each sample is closest to
        for sample in list(self.ctDf):
            temp = 0
            cEle = ''
            for cType in list(self.cDf):    
                if abs(self.cDf[cType].corr(self.ctDf[sample], method = method)) > temp:
                    temp = abs(self.cDf[cType].corr(self.ctDf[sample], method = method))
                    cEle = cType
            cList.append(cEle)
        #convert the list of cancer types to colors, then joining it with the col_colors
        cDict = {'bronchioid':'Red', 'magnoid':'Blue' , 'squamoid':'Yellow'}
        centroidClusters = [inner for outer in [[counter+1 for i in x] for counter, x in enumerate([x  for x in self.sortedClusters]) ] for inner in outer]
        
        c2List = [cDict[x] for x in cList]
        self.centroid = pd.DataFrame(c2List, columns = ['Centroid'])
        self.centroid.index = self.tDf.index
        self.centroid['Clusters'] = centroidClusters
        self.centroid['Cancers'] = cList
        self.col_colors = pd.concat([self.col_colors, self.centroid['Centroid']], axis = 1)
        ax = sns.countplot(x="Clusters", hue="Cancers", data=self.centroid)
                
        
    def logFoldChange(self, csvpath = 'data/test.csv'):
        self.readCsv(csvpath = csvpath, readType = 'lfc')
        self.lfcDf = pd.DataFrame()
        for counter, columns in enumerate(self.sortedClusters):
            self.lfcDf = pd.concat([self.lfcDf, self.rawDf.filter(items = columns).mean(axis = 1)],axis  =1)
        self.lfcDf.columns = [f'Cluster {cluster}' for cluster in range(1, 1+len(self.sortedClusters) )]
        tempDf = self.lfcDf.copy()
        for i, column in enumerate(list(self.lfcDf)):
            self.lfcDf[f'Cluster {i+1}'] = np.log2(tempDf.loc[:, tempDf.columns != f'Cluster {i+1}'].mean(axis=1).div(tempDf[f'Cluster {i+1}']).replace(0,1)).fillna(0)
            self.lfcDf[f'Cluster {i+1}'] = np.log2(tempDf[f'Cluster {i+1}'].div(tempDf.loc[:, tempDf.columns != f'Cluster {i+1}'].mean(axis=1)).replace(0,1)).fillna(0)
        newIndex = [x.split('.')[0] for x in list(self.lfcDf.index)]
        self.lfcDf.index = newIndex
        
        
    def normalityTest(self):
        pvals = list()
        for column in tqdm(list(self.tDf)):
            pvals.append(normaltest(self.tDf[column])[1])
        sns.displot(pvals, binwidth = 0.05)
        
    def kaplanMeier(self):
        if self.clinicalDf is None:
             self.readCsv(readType = 'Clinical')
        df = self.clinicalDf
        km = kmf() 
        
        for counter, cluster in enumerate(self.sortedClusters):
            patients = ['-'.join(x.split('-')[0:3]) for x in cluster]
            columns = ['vital_status', 'days_to_death', 'days_to_last_followup']
            df.filter(items = patients, axis = 'index')[columns]
            f = lambda x: 1 if x == 'Alive' else  0
            status = [f(x) for x in df.filter(items = patients, axis = 'index')['vital_status'].tolist()]
            followup = df.filter(items = patients, axis = 'index')['days_to_last_followup'].tolist()
            status = df.filter(items = patients, axis = 'index')['vital_status'].tolist()
            death = df.filter(items = patients, axis = 'index')['days_to_death'].tolist()
            Died, time = list(), list()
            for i, stat in enumerate(status):
                if stat == 'Alive':
                    Died.append(1)
                    time.append(followup[i])
                else:
                    try:
                        float(death[i])
                        Died.append(0)
                        time.append(int(float(death[i])))
                    except:
                        continue
            time = [(float(x)/365) for x in time]
            km.fit(time,Died,  label=f'Cluster {counter+1}')
            if counter == 0:
                ax = km.plot(ci_show=False, show_censors = True)
            elif counter % 2 != 0:
                ax1 = km.plot(ax=ax, ci_show=False, show_censors = True)
            else:
                ax = km.plot(ax=ax1, ci_show=False, show_censors = True)  
    def coxHazard(self):
        if self.clinicalDf is None:
             self.readCsv(readType = 'Clinical')
        
        
    def addFeature(self, feature = None, attr = None, censor = None, title = None):
        if title is None: title = feature
        if self.clinicalDf is None:
             self.readCsv(readType = 'Clinical')
        patients = ['-'.join(x.split('-')[0:3]) for x in self.allClusters]
        df = self.clinicalDf.filter(items = patients, axis = 'index')
        df.index = [f'{x}-01A' for x in df.index]
        if censor is None: s = df[feature].replace(regex={r'^((?!{attr}).*)'.format(attr = attr): 'silver', r'^{attr}'.format(attr = attr): 'black'})
        else: 
            if censor[0] == '[':
                censor = '\\'+censor
            s = df[feature].replace(regex={r'^((?!{attr}|{cens}).*)'.format(attr = attr, cens = censor): 'silver', r'^{attr}'.format(attr = attr): 'black', r'^{cens}'.format(cens = censor): 'white' })
        s = s.rename(title)
        self.col_colors = pd.concat([self.col_colors, s], axis = 1)
    
    def dropFeature(self, feature):
        self.col_colors.drop(columns= ['Non-smokers'])
        
    def plotSingleGene(self, genes = list()):
        if self.symbolDf is None: self.ensembl2symbol()
        fig, axs = plt.subplots(len(genes), figsize = (10,10))
        for counter, gene in enumerate(genes):
            lista = [gene, 'Clusters']
            g  = sns.violinplot(data = self.tSymbolDf[lista], x='Clusters', y =gene, ax = axs[counter])
            g.axhline(0)
        fig.tight_layout()
            
    def plotMultipleGenes(self, geneList = list):
        if self.symbolDf is None: self.ensembl2symbol()
        if len(geneList) > 1: fig, axs = plt.subplots(len(geneList), figsize = (10,10))
        for counter, genes in enumerate(geneList):
            t = pd.melt(self.tSymbolDf.filter(items = genes+ ['Clusters']), id_vars = 'Clusters' ).drop(['variable'], axis = 1)
            if len(geneList) > 1: g =  sns.boxplot(data = t, x = 'Clusters', y = 'value', ax = axs[counter])
            else: g = sns.violinplot(data = t, x = 'Clusters', y = 'value')
            g.axhline(0)
        if len(geneList) > 1: fig.tight_layout()
            
    def readData(self):
        samples = ['-'.join(sample.split('-')[0:3]) for sample in self.df.columns.tolist()]
        sortedL = [['-'.join(x.split('-')[0:3]) for x in cluster] for cluster in self.sortedClusters]
        f = lambda l: ['-'.join(x.split('-')[0:3]) for x in l]
        fDot = lambda l: ['-'.join(x.split('.')[0:3]) for x in l]
        
        
        neoDf = pd.read_csv('data/clinicalData/TCGA_PCA.mc3.v0.2.8.CONTROLLED.filtered.sample_neoantigens_10062017.tsv', sep = '\t')
        neoDf = neoDf[neoDf['sample'].isin(samples)]
        neoDf = neoDf.set_index('sample')
        
        pmhcDf = pd.read_csv('data/clinicalData/TCGA_pMHC_SNV_sampleSummary_MC3_v0.2.8.CONTROLLED_170404.tsv', '\t')
        pmhcDf['barcode'] = f(pmhcDf['barcode'].tolist())
        pmhcDf = pmhcDf[pmhcDf['barcode'].isin(samples)]
        pmhcDf = pmhcDf.set_index('barcode')
        
        mutDf = pd.read_csv('data/clinicalData/mutation-load_updated.txt', sep='\t')
        mutDf = mutDf[mutDf['Patient_ID'].isin(samples)]
        mutDf = mutDf.set_index('Patient_ID')
        
        
        mastDf = pd.read_csv('data/clinicalData/TCGA_mastercalls.abs_tables_JSedit.fixed.txt',sep = '\t')
        mastDf['sample'] = f(mastDf['sample'].tolist())
        mastDf = mastDf[mastDf['sample'].isin(samples)]
        mastDf = mastDf.set_index('sample')
        
        absDf = pd.read_csv('data/clinicalData/ABSOLUTE_scores.tsv', sep='\t')
        absDf['Unnamed: 0'] = f(absDf['Unnamed: 0'].tolist())
        absDf = absDf[absDf['Unnamed: 0'].isin(samples)]
        absDf = absDf.set_index('Unnamed: 0')
        
        leuDf = pd.read_csv('data/clinicalData/TCGA_all_leuk_estimate.masked.20170107.tsv', sep ='\t')
        leuDf.columns = ['Type', 'Sample', 'Leukocyte Fraction']
        leuDf['Sample'] = f(leuDf['Sample'].tolist())
        leuDf = leuDf[leuDf['Sample'].isin(samples)]
        leuDf = leuDf.set_index('Sample')
        leuDf = leuDf[~leuDf.index.duplicated(keep='first')]
        
        cibDf = pd.read_csv('data/clinicalData/TCGA.Kallisto.fullIDs.cibersort.relative.tsv', sep = '\t')
        cibDf['SampleID'] = fDot(cibDf['SampleID'].tolist())
        cibDf = cibDf[cibDf['SampleID'].isin(samples)]
        cibDf = cibDf[~cibDf['SampleID'].duplicated(keep='first')]
        cibDf = cibDf.set_index('SampleID')
        
        mitDf = pd.read_csv('data/clinicalData/mitcr_sampleStatistics_20160714.tsv', sep='\t')
        mitDf = mitDf[mitDf['ParticipantBarcode'].isin(samples)]
        mitDf = mitDf.set_index('ParticipantBarcode')
        mitDf = mitDf[~mitDf.index.duplicated(keep='first')]
        
        subDf = pd.read_csv('data/clinicalData/TCGASubtype.20170308.tsv', sep ='\t')
        subDf = subDf[subDf['pan.samplesID'].isin(samples)]
        subDf = subDf.set_index('pan.samplesID')
        
        hrdDf = pd.read_csv('data/clinicalData/TCGA.HRD_withSampleID.txt',sep = '\t')
        hrdDf['sampleID'] = f(hrdDf['sampleID'].tolist())
        hrdDf = hrdDf[hrdDf['sampleID'].isin(samples)].set_index('sampleID')
        
        self.eDf = pd.concat([leuDf, cibDf, mastDf, absDf, pmhcDf, neoDf, hrdDf, subDf], axis = 1)
        temp = self.tDf['Clusters']
        temp.index = f(temp.index.tolist())
        self.eDf['Clusters'] = temp
        
    def corrDNA(self):
        if self.eDf is None: self.readData()
        f = lambda l: ['-'.join(s.split('-')[0:3]) for s in l]
        hDf = pd.DataFrame()
        for i, samples in enumerate(self.sortedClusters):
            hDf[f'Cluster {i+1}'] = self.eDf.filter(items=f(samples), axis = 'index').corr()['Leukocyte Fraction']
        hDf = hDf.dropna()
        self.hDf = hDf.iloc[1:]
        g = sns.clustermap(self.hDf, cmap = 'vlag', yticklabels = True, row_cluster = False)

        
    def plotSignatures(self):
        df = pd.read_csv('data/clinicalData/signature_profile_sample_mSignatureDB.txt', sep = '\t')
        df['Tumor_Sample_Barcode'] = [x.replace('.', '-') for x in df['Tumor_Sample_Barcode'].tolist()]
        df = df[df['Tumor_Sample_Barcode'].isin(['-'.join(sample.split('-')[0:3]) for sample in self.samples])]
        newC = self.col_colors
        newC.index = ['-'.join(x.split('-')[0:3]) for x in list(newC.index)]
        self.sigDf = df.pivot(index='Signature', columns = 'Tumor_Sample_Barcode', values = 'Contribution')[self.samples]
        sns.clustermap(self.sigDf, col_colors = newC, row_cluster = True, col_cluster = False)
    def boxplotExternalData(self, cols = 6):
        if self.eDf is None: self.readData()
        df = self.eDf
        for column in list(df):
            if not is_numeric_dtype(df[column]): df = df.drop([column], axis = 1)
        df['id'] = df.index.tolist()
        self.newDf = pd.melt(df, id_vars = ['id', 'Clusters'])
        g = sns.catplot(data = self.newDf, kind = 'box', col = 'variable', y ='value', x = 'Clusters', col_wrap = cols, sharey = False)