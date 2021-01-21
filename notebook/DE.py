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
  
    def readCsv(self, csvpath, readType = 'csv', newIndex = False, sep = ',', fpkm = False):
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
                self.metaDf = pd.read_excel(r'130932_1_supp_2656123_nbr7l8.xlsx')
                self.metaDf = self.metaDf.iloc[:, 0:2]
                
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
        #Reads the output of SRIQ clustring
        self.df = pd.read_csv(csvpath, delimiter = '\t')
        self.df = self.df.set_index([columnname])
        index = list(self.df.index)
        index.append('Clusters')
        self.clusters = dict()
        self.allClusters = list()
        #Gathers all samples in allClusters and specific clusters in clusters
        for filename in os.listdir(clusterpath):
            self.key = filename[1:]
            string = ''
            with open(os.path.join(clusterpath, filename)) as file:
                for line in file.readlines():
                    string = string + line
                    if line.strip() != columnname: self.allClusters.append(line.strip())
                string.replace(columnname, '')
                string = string.split()
                string.pop(0)
                self.clusters[filename] = string
                self.allClusters = list(dict.fromkeys(self.allClusters))
        self.df = self.df.filter(items=self.allClusters)
        dicts = dict()
        
        #Creates the clusters row to be appended into the working dataframe
        for i in range(len(self.clusters)):
            for sample in self.clusters['{}{}'.format(i+1, self.key)]:
                dicts[sample] = i+1
        
        self.df = self.df.append(dicts, ignore_index= True)
        
        self.df = self.df.set_index([index])
        
        #Sorts the columns so that cluster 1 is first and so forth
        newDf = pd.DataFrame()
        for i in range(len(self.clusters)):
            for sample in self.clusters['{}{}'.format(i+1, self.key)]:
                newDf[sample] = self.df[sample]
        newDf = newDf.set_index([index])
        self.df = newDf 
        self.tDf = self.df.transpose()
        #Creates a dataframe with the colors of the clusters to be added for the clustermap
        lut = dict(zip(self.tDf['Clusters'].unique(), 'bgrcmykw'))
        self.col_colors = self.tDf['Clusters'].map(lut)
        ax = sns.countplot(self.tDf['Clusters'])
        self.df = self.df[(self.df.T != 0).any()]
        self.tDf = self.df.transpose()
        self.filterDf = self.df
        self.sortedClusters = [self.clusters[f'{x+1}_{clusterpath[-1]}C_Samples.txt'] for x in range(len(self.clusters))]
        
        
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
        for i in range(len(self.clusters)):
            self.pDf['Cluster {} corr'.format(i+1)] = sm.stats.multitest.fdrcorrection(self.pDf['Cluster {}'.format(i+1)].array, alpha=0.05, method='indep', is_sorted=False)[1]
        print('Filtering out genes')
        self.pDf = self.pDf.loc[(self.pDf.iloc[:,4:] < 0.01).sum(axis = 1) <4]
        bla = pd.DataFrame(columns = [f'Cluster {i+1} corr:' for i in range(len(self.sortedClusters))])
        # Sets the lesser significant genes to a p value of 1
        tempDf = self.pDf.iloc[:,len(self.clusters):].min(axis = 1)
        for i in tqdm(range(self.pDf.shape[0])):
            bla = bla.append(self.pDf.iloc[:,len(self.clusters):].iloc[i].transform(lambda x: x if x == tempDf.iloc[i] else x/x))
        self.pDf.iloc[:,len(self.clusters):] = bla
        
            
    def filtering(self, threshold = 2, filteringType = True, start = 1 ,transformed = True, csvpath = '//home/jacob/examensarbete/test.csv', sep = '.'):
        #Filters based on fpkm value. Not needed if using the test metho use threshold = 1
        if self.pDf is not None:
            if filteringType == 'fpkm': self.readCsv(csvpath, readType = 'fpkm',  sep = ',', fpkm = True)
            elif filteringType == 'log2fold': self.logFoldChange()
            #
            genes = list()
            #Fethes the significant genes
            for i in range(len(self.clusters)):
                genes.append(self.pDf[self.pDf['Cluster {} corr'.format(i+1)] < 0.01].index.tolist())
            self.tList = list()
            #Filters the genes based of fpkm
            if filteringType == 'fpkm':
                for i in range(len(self.clusters)):
                    finalList = genes[i]
                    lista = list()
                    for gene in finalList:
                        if  (self.fpkmDf[gene] > threshold):
                            lista.append(gene)
                    self.tList.append(lista)
            if filteringType == 'log2fold':
                for i in range(len(self.clusters)):
                    finalList = genes[i]
                    #I dont fucking know
                    self.tList.append(list(set(genes[i]) & set(list(self.lfcDf.loc[abs(self.lfcDf[f'Cluster {i+1}']) > threshold].index))))
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
        
    
    def enrichR(self, cluster, dbs):
        #Fetches enrichment info from enrichR against different databases.
        #Fetches user list id for api calls
        if self.eList is not None:
            ENRICHR_URL = 'http://maayanlab.cloud/Enrichr/addList'
            genes_str = '\n'.join(self.eList[cluster])
            description = 'Example gene list'
            payload = {
                'list': (None, genes_str),
                'description': (None, description)
            }
            response = requests.post(ENRICHR_URL, files=payload)
            if not response.ok:
                raise Exception('Error analyzing gene list')
            user_list_id = json.loads(response.text)['userListId']
            c = [0,1,2,3,4,5,6,7,8,'db']
            self.rDf = pd.DataFrame(columns = c)
            #Fetches values from each database and saved them longformat into rDf
            for db in dbs:
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
                testDf['db'] = [db]*5
                self.rDf = self.rDf.append(testDf)
                self.rDf = self.rDf.sort_values(by = [2], ascending = False)
        else: print('Error: Run the filtering module')
            
    def plotBar(self):
        #plots the enirchment result from rDf
        if self.rDf is not None:
            g = sns.catplot(x= 2, y=1,kind='bar', data = self.rDf)
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
        g = sns.catplot(x= 'num', y='per_var',kind='bar', data = self.sDf, col='cutoff',col_wrap= 2)
        g2 = sns.relplot(x= 'comp 1', y='comp 2', data= self.uDf, style = 'hue', hue='hue', col='cutoff',col_wrap=2,facet_kws=dict(sharex=False))
        
        #Adjusts the axices the plots
        for i, ele in enumerate(lista):
            print(g2.axes[i])
            temp = self.uDf.loc[self.uDf['cutoff'] == ele]
            xmin, xmax, ymin, ymax = temp['comp 1'].min(), temp['comp 1'].max(), temp['comp 2'].min(), temp['comp 2'].max()
            g2.axes[i].set_xlim(xmin-1, xmax+1)
            g2.axes[i].set_ylim(ymin-1, ymax+1)
            print(g2.axes[i])
        print('done')
        
    def ensemble2gene(self):
        #Converts the enrichment list from ensembl to gene symbol
        mg = mygene.MyGeneInfo()
        for i, ele in enumerate(self.eList):
            templist = list()
            listan = mg.querymany(ele, scopes = 'ensembl.gene')
            for dicts in listan:
                if 'symbol' in dicts:
                    templist.append(dicts['symbol'])
            self.eList[i] = templist    

    def metaGenes(self):
        df = pd.read_excel(r'130932_1_supp_2656123_nbr7l8.xlsx')
        df = df.iloc[:, 0:2]
        boxDf = pd.DataFrame(columns = ['Values', 'Cluster', 'Metagene'])
        for x in df['Network cluster'].unique():
            print(x)
            for i, cluster in enumerate(self.sortedClusters):
                values = (self.symbolDf[cluster].filter(items = list(df[df['Network cluster'] == x]['Gene Symbol']), axis = 0).mean(axis = 1).values)
                tempDf = pd.DataFrame({'Values': values, 'Cluster':[i+1 for c in range(len(values))], 'Metagene':[x for c in range(len(values))]})
                boxDf = pd.concat([boxDf, tempDf])
        ax = sns.catplot(kind = 'box', data = boxDf, x = 'Cluster', y = 'Values', col = 'Metagene', col_wrap = 2)
    
        
        
        
    def ensembl2symbol(self):
        print('This may produce duplicate indexes')
        newIndex = [x.split('.')[0] for x in list(self.df.index)]
        self.df.index = newIndex
        self.tDf = self.df.transpose()
        mg = mygene.MyGeneInfo()
        print('Finding symbols for ensemblIds')
        genes = mg.querymany(newIndex, scopes = 'ensembl.gene', returnall= True, as_dataframe= True, verbose = False)
        genes = genes['out']
        newIndex = (genes['symbol'].dropna())
        self.symbolDf = self.df.reindex(list(newIndex.index))
        self.symbolDf.index = list(newIndex)
        self.symbolDf = self.symbolDf[~self.symbolDf.index.duplicated(keep='first')]

        
                
        
    def centroids(self, centroidpath = '/home/jacob/examensarbete/wilkerson.2012.LAD.predictor.centroids.csv', ensembl = True, method = 'pearson'):
        #Purpose: calculating centroids on cancer data
        #read the centroid predictor data set
        self.readCsv(csvpath = centroidpath, readType = 'centroid')
        if ensembl:
            #If the working data is in ensemble, it will find the matchning gene-ensembldids.
            #Finallist consist of the matching genes
            if self.symbolDf is None:
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
        centroidClusters = [inner for outer in [[counter+1 for i in x] for counter, x in enumerate([self.clusters[f'{x}_4C_Samples.txt'] for x in range(1,1+len(self.clusters))]) ] for inner in outer]
        
        c2List = [cDict[x] for x in cList]
        self.centroids = pd.DataFrame(c2List, columns = ['Colors'])
        self.centroids.index = self.tDf.index
        self.centroids['Clusters'] = centroidClusters
        self.centroids['Cancers'] = cList
        self.col_colors = pd.concat([self.col_colors, self.centroids['Colors']], axis = 1)
        ax = sns.countplot(x="Clusters", hue="Cancers", data=self.centroids)
                
        
    def logFoldChange(self, csvpath = '//home/jacob/examensarbete/test.csv'):
        self.readCsv(csvpath = csvpath, readType = 'lfc')
        self.lfcDf = pd.DataFrame()
        for counter, columns in enumerate(self.sortedClusters):
            self.lfcDf = pd.concat([self.lfcDf, self.rawDf.filter(items = columns).mean(axis = 1)],axis  =1)
        self.lfcDf.columns = [f'Cluster {cluster}' for cluster in range(1, 1+len(self.clusters) )]
        tempDf = self.lfcDf.copy()
        for i, column in enumerate(list(self.lfcDf)):
            self.lfcDf[f'Cluster {i+1}'] = np.log2(tempDf.loc[:, tempDf.columns != f'Cluster {i+1}'].mean(axis=1).div(tempDf[f'Cluster {i+1}']).replace(0,1)).fillna(0)
        newIndex = [x.split('.')[0] for x in list(self.lfcDf.index)]
        self.lfcDf.index = newIndex
        
        
    def normalityTest(self):
        pvals = list()
        for column in tqdm(list(self.tDf)):
            pvals.append(normaltest(self.tDf[column])[1])
        sns.displot(pvals, binwidth = 0.05)
        
    def kaplanMeier(self):
        df = pd.read_csv(r'clinical_PANCAN_patient_with_followup.tsv', sep = '\t', engine = 'python', index_col = 'bcr_patient_barcode')

        km = kmf() 
        
        for counter, cluster in enumerate(self.sortedClusters):
            patients = ['-'.join(x.split('-')[0:3]) for x in cluster]
            columns = ['vital_status', 'days_to_death', 'days_to_last_followup']
            df.filter(items = patients, axis = 'index')[columns]
            # def f(string):
            #     if string == 'Alive': return 1
            #     else: return 0
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
                ax = km.plot(ci_show=False)
            elif counter % 2 != 0:
                ax1 = km.plot(ax=ax, ci_show=False)
            else:
                ax = km.plot(ax=ax1, ci_show=False)  
    def coxHazard(self):
        df = pd.read_csv(r'clinical_PANCAN_patient_with_followup.tsv', sep = '\t', engine = 'python', index_col = 'bcr_patient_barcode')
        df = df[self.sortedClusters[0]]
        
        
        
