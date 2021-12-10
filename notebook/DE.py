#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 09:04:10 2020

@author: jacob karlstrom
"""
#import openpyxl
import json
import requests
import seaborn as sns; sns.set_theme()
import os
import sys
import math
from numpy import *
import sklearn
from sklearn.metrics import silhouette_samples, silhouette_score
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd
import seaborn as sns
import umap.umap_ as umap
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
from matplotlib.colors import ListedColormap



class networkAnalysis():
    def __init__(self, newDf = None, symbol = False):
        #Initiates the working dataframe for manual initiation, not recomended.
        if newDf:
            self.gexDf = newDf
            self.transposedGexDf = self.gexDf.transpose()
        #Initialized the different dataframes used in the object
        self.eList = None
        self.pValuesDf = None
        self.lfcDf = None
        self.goEnrichDf = None
        self.symbolDf = None
        self.clinicalDf = None
        self.metaGeneDf = None
        self.immuneDf = None
        self.pca_data = None
        self.samDf = None
        self.preDf = None
        self.samResults = None
        
    def SilhouttePlot(self, sizes = 100):
        #if self.pca_data is None:
        #    self.calcPcaModel()
        if self.samResults is not None:
            self.calcPcaModel()
            labels = self.gexDf.iloc[-1:].transpose()['Clusters']
            X = self.samResults.corr('pearson')
            X = X.applymap(lambda x: (1 - x))
        else:
            self.calcPcaModel()
            labels = self.gexDf.loc['Clusters'] 
            X = self.gexDf.drop('Clusters', axis = 0).corr('pearson')
            X = X.applymap(lambda x: (1 - x))
        
        np.fill_diagonal(X.to_numpy(), 0)
        
        
        score = sklearn.metrics.silhouette_score(X, labels, metric='precomputed')
        
        samples = silhouette_samples(X, labels)
        
        fig, (ax1, ax2) = plt.subplots(1, 2)
        fig.set_size_inches(18, 7)
        y_lower = 10
        colors = sns.color_palette(n_colors = self.clusterNum)
        cs = sns.color_palette(as_cmap = True)
        my_cmap = ListedColormap(sns.color_palette(cs).as_hex())
        colors2 = [[cs[i] for x in range(len(self.sortedClusterList[i]))] for i in range(self.clusterNum)]
        colors3 = list()
        [[colors3.append(x) for x in l] for l in colors2]
        for i in range(1, int(self.clusterNum)+1):
            # Aggregate the silhouette scores for samples belonging to
            # cluster i, and sort them
            ith_cluster_silhouette_values = \
                samples[labels == i]
        
            ith_cluster_silhouette_values.sort()
        
            size_cluster_i = ith_cluster_silhouette_values.shape[0]
            y_upper = y_lower + size_cluster_i
            color = colors[i-1]
            ax1.fill_betweenx(np.arange(y_lower, y_upper),
                                0, ith_cluster_silhouette_values,
                                facecolor=color, edgecolor=color, alpha=0.7)
        
            # Label the silhouette plots with their cluster numbers at the middle
            ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))
        
            # Compute the new y_lower for next plot
            y_lower = y_upper + 10  # 10 for the 0 samples
            
        ax1.set_title("The silhouette plot for the various clusters.")
        ax1.set_xlabel("The silhouette coefficient values")
        ax1.set_ylabel("Cluster label")
        #self.colors2 = cm.nipy_spectral(labels.astype(float) / self.clusterNum)
        ax2.scatter(self.pca_data[:, 0], self.pca_data[:, 1], marker='.', s=sizes, lw=0, alpha=0.7,
                    c=colors3,cmap = my_cmap, edgecolor='k')
        #centers = clusterer.cluster_centers_
        
        ax2.set_title("The visualization of the clustered data.")
        ax2.set_xlabel("Feature space for the 1st feature")
        ax2.set_ylabel("Feature space for the 2nd feature")
        
        # The vertical line for average silhouette score of all the values
        ax1.axvline(x=score, color="red", linestyle="--")
        
        ax1.set_yticks([])  # Clear the yaxis labels / ticks
        ax1.set_xticks([x for x in np.arange(round(min(samples), 1), round(max(samples), 1), 0.1)])
        plt.suptitle(("Silhouette analysis for SRIQ clustering on sample data "
                  "with n_clusters = %d" % self.clusterNum),
                 fontsize=14, fontweight='bold')
        plt.show()
        
    def preProcess(self):

        DEGDf = self.gexDf
        DEGDf = DEGDf[(DEGDf.iloc[:,1:].T != 0).any()]
        DEGDf.iloc[:,1:] = DEGDf.iloc[:,1:].apply(lambda x: x+0.1)
        DEGDf.iloc[:,1:] = DEGDf.iloc[:,1:].apply(lambda x: x/x.median())
        DEGDf.iloc[:,1:] = np.log2(DEGDf.iloc[:,1:])
        self.DEGDf = DEGDf


        f = lambda x: 1 if x < 1 else x
        self.gexDf.iloc[:,1:] = self.gexDf.iloc[:,1:].applymap(func = f)
        self.preDf = self.gexDf
        #self.gexDf.iloc[:,1:] = self.gexDf.iloc[:,1:].applymap(lambda x: x + 0.01)
        self.gexDf.iloc[:,1:] = self.gexDf.transpose().iloc[1:].transform(lambda x: x/x.median()).transpose()
        self.preFilterCalc()
        self.gexDf.iloc[:,1:] = self.gexDf.iloc[:,1:].applymap(lambda x: np.log2(x))
        self.gexDf = self.gexDf[(self.gexDf.iloc[:,1:].T != 0).any()]
        
        self.DEGDf =  self.DEGDf[self.DEGDf['Gene'].isin(self.gexDf['Gene'].tolist())]

        
    def preFilterCalc(self):
        self.fpkmVarDf = self.gexDf.iloc[:,1:].var(axis = 'columns')
        self.fpkmVarDf.index = self.gexDf['Gene'].tolist()
        self.fpkmVarDf = self.fpkmVarDf.sort_values(ascending = True)
        
    def preFilter(self, bottom):    
        self.filterDf = self.gexDf[self.gexDf['Gene'].isin(self.fpkmVarDf.iloc[int(len(self.fpkmVarDf)*bottom):].index.tolist())]
        
    def DEGFile(self, filename = None):
        df = self.gexDf
        f = lambda x: 1 if x < 1 else x
        df.iloc[:,1:] = df.iloc[:,1:].applymap(func = f)
        df = np.log2(df)
        df = df[(df.T != 0).any()]
        if filename is None: df.to_csv(f'data/DEG.csv')
    
  
    def readCsv(self, csvpath = None, readType = 'csv', newIndex = False, sep = ',', fpkm = False):
            if readType == 'csv':
                self.gexDf = pd.read_csv(csvpath, sep = sep, index_col=False)
                self.transposedGexDf = self.gexDf.transpose()  
            elif readType == 'fpkm':
                self.fpkmDf = pd.read_csv(csvpath, sep = sep).set_index('Unnamed: 0').transpose()
                meanList = list()
                for gene in list(self.fpkmDf):
                    meanList.append(self.fpkmDf[gene].mean())
                self.fpkmDf = self.fpkmDf.transpose()
                self.fpkmDf['FPKM'] = meanList
                self.fpkmDf = self.fpkmDf['FPKM']
            elif readType == 'centroid':
                self.cDf = pd.read_csv(csvpath, sep = sep, index_col = 'Unnamed: 0',error_bad_lines=False)
            elif readType=='lfc':
                self.rawDf = pd.read_csv(csvpath, sep = sep)
            elif readType == 'meta-gene':
                self.metaGeneDf = pd.read_csv(r'data/extraData/metagenes2.csv', error_bad_lines=False, sep = ';')
                self.metaGeneDf = self.metaGeneDf.iloc[:, 0:2]
            elif readType == 'Clinical':
                if csvpath is None:
                    self.clinicalDf = pd.read_csv(r'data/extraData/clinical_PANCAN_patient_with_followup.tsv', sep = '\t', engine = 'python', index_col = 'bcr_patient_barcode')
                else:
                    self.clinicalDf = pd.read_csv(csvpath, sep = '\t', engine = 'python', index_col = 'bcr_patient_barcode')
            if newIndex:
                self.gexDf = self.gexDf.set_index(newIndex)
                self.transposedGexDf = self.gexDf.transpose()
                
            
    def filterCCL(self, SRIQpath):
        lista = list()
        for file in sorted(os.listdir(SRIQpath)):
            lista.append([line.strip() for line in open(f'{SRIQpath}/{file}')][1:])
        lista2 = list()
        [[lista2.append(x) for x in l] for l in lista]
        self.gexDf = self.gexDf.filter(items = lista2, axis = 'columns')
        self.transposedGexDf = self.gexDf.transpose()
        self.sortedClusterList = [[x for x in l if x in lista2] for l in self.sortedClusers]

        
        
    def filterVariantGenes(self, start = 0, top = 1, bottom = 0):
        #Filter the least variant genes based on threshold ranging 0 to 1
        self.filterDf = self.gexDf.iloc[:,start:].loc[(self.gexDf.iloc[:,start:].var(axis = 1) < self.gexDf.iloc[:,start:].var(axis = 1).quantile(top)) & (self.gexDf.iloc[:,start:].var(axis = 1) > self.gexDf.iloc[:,start:].var(axis = 1).quantile(bottom))]
    
    def readSRIQ(self, csvpath, clusterpath, columnname = 'Gene', symbol = False, **kwargs):
        self.gexDf = pd.read_csv(csvpath, delimiter = '\t', index_col = columnname)
        self.sortedClusterList = list()
        self.clusterList = list()
        for file in sorted(os.listdir(clusterpath)):
            lista = [line.strip() for line in open(f'{clusterpath}/{file}')][1:]
            self.clusterList += lista
            self.sortedClusterList.append(lista)
        self.gexDf = self.gexDf[self.clusterList]
        clusterRow = [inner for outer in [[counter+1]*len(x) for counter, x in enumerate(self.sortedClusterList)] for inner in outer]
        clusterRow = pd.Series(data = clusterRow, name = 'Clusters', index = self.gexDf.columns)
        self.gexDf = self.gexDf.append(clusterRow)
        self.gexDf = self.gexDf[(self.gexDf.T != 0).any()]
        self.transposedGexDf = self.gexDf.transpose()
        lut = dict(zip(self.transposedGexDf['Clusters'].unique(), sns.color_palette(n_colors = len(os.listdir(clusterpath)),as_cmap = True)))
        self.col_colors = self.transposedGexDf['Clusters'].map(lut)
        ax = sns.countplot(self.transposedGexDf['Clusters'], **kwargs)
        self.filterDf = self.gexDf
        self.clusterNum = len(os.listdir(clusterpath))
        f = lambda l: ['-'.join(s.split('-')[0:3]) for s in l]
        self.samples = f(self.gexDf.columns.tolist())
            
        self.dist = clusterpath
        if symbol:
            self.symbolDf = self.gexDf
            self.tSymbolDf = self.symbolDf.transpose()
    def setlog2Df(self, df):
        self.gexDf = df
        self.transposedGexDf = self.gexDf.transpose()
        

    def calcPcaModel(self, Clustered = False, demo = False):
        #Calculates the pca model and plots the scree plot
        pca = PCA()
        #if Clustered: 
            #pca.fit(self.transposedGexDf.drop(['Clusters'],axis=1))
            #self.pca_data = pca.transform(self.transposedGexDf.drop(['Clusters'],axis=1))
        if self.samResults is not None:
            pca.fit(self.samResults.transpose())
            self.pca_data = pca.transform(self.samResults.transpose())
        else: 
            pca.fit(self.transposedGexDf.drop(['Clusters'],axis=1))
            self.pca_data = pca.transform(self.transposedGexDf.drop(['Clusters'],axis=1))
            #pca.fit(self.transposedGexDf)
            #self.pca_data = pca.transform(self.transposedGexDf)
        
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
        


        
    def pcaPlot(self, Clustered = True, **kwargs):
        #Scatterplots the pca model
        self.pca_df = pd.DataFrame(self.pca_data,  columns=self.labels)
        if Clustered:
            self.pca_df = self.pca_df.join(self.transposedGexDf.filter(['Clusters']).reset_index())
            sns.relplot(x="PC1", y="PC2", data=self.pca_df, hue='Clusters', style='Clusters', **kwargs);
        else:
            sns.relplot(x="PC1", y="PC2", data=self.pca_df, **kwargs);
        plt.show()
        
        
    def Umap(self,n_neigh = 10, min_dist = 0.3, n_comp = 4, demo = False, **kwargs):
        #Does a umap model based on the pca model, n = number of desired principal components
        reducer = umap.UMAP(n_neighbors=n_neigh, min_dist= min_dist,  n_components = n_comp)

        embedding = reducer.fit_transform(self.gexDf.iloc[:-1,:].transpose())
        

        sns.scatterplot(x = embedding[:,0], y = embedding[:,1], hue = self.gexDf.loc['Clusters'], palette = sns.color_palette(n_colors = self.clusterNum))

        plt.xlabel('UMAP1')
        plt.ylabel('UMAP2')

        plt.title(f'UMAP of SRIQ K{self.clusterNum} solution')

    def configResources(self,outPath = '/Users/jacobkarlstrom/projekt/SRIQ/software/output/',studyPath='/Users/jacobkarlstrom/projekt/SRIQ/notebook/data/expressionData/',  data='filtered(21k)', resources = '../software/VRLA/resources/test.properties', studyName = 'SRIQ', cutOff = None, permutations = 10000, iterations = 10, minBagSize=1200, minClusterSize = 0):
        if cutOff is None: cutOff = [0.9,0.89,0.88,0.87,0.86,0.85,0.84,0.83,0.82,0.81,0.80,0.79,0.78,0.77,0.76,0.75,0.74,0.73,0.72,0.71,0.7,0.69,0.68,0.67,0.66,0.65,0.64,0.63,0.62,0.61,0.6,0.59,0.58,0.57,0.56,0.55,0.54,0.53,0.52,0.51,0.5,0.49,0.48,0.47,0.46,0.45,0.44,0.43,0.42,0.41,0.4,0.39,0.38,0.37,0.36,0.35,0.34]
        output = ''
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
                elif 'studyPath' in line:
                    output += f'studyPath={studyPath}\n'
                elif 'outPath' in line:
                    output += f'outPath={outPath}\n'
                else:
                    output += line
        with open(resources, 'w') as file:
            file.write(output)
    def diffGeneAnalysis(self, transformed = True, test = 'non-parametric'):
        
        #Tests each gene for its significance in the different clusters
        self.pValuesDf = pd.DataFrame(index = [x.split('.')[0] for x in self.filterDf.index.tolist()])
        if test == 't-test':
            for counter, columns in enumerate(tqdm(self.sortedClusterList)):
                res = list()
                set2 = set(self.filterDf) - set(columns)
                df1, df2 = self.filterDf[columns], self.filterDf[set2]
                for i in tqdm(range(df1.shape[0])):
                    res.append(ttest_ind(df1.iloc[i].values, df2.iloc[i].values,equal_var = False)[1])
                    if np.isnan(res[-1]):
                        res.pop(-1)
                        res.append(1)
                self.pValuesDf[f'Cluster {counter +1}'] = res
        elif test == 'mannwhitneyu':
            for counter, columns in enumerate(tqdm(self.sortedClusterList)):
                res = list()
                set2 = set(self.filterDf) - set(columns)
                df1, df2 = self.filterDf[columns], self.filterDf[set2]
                for i in tqdm(range(df1.shape[0])):
                    res.append(mannwhitneyu(df1.iloc[i].values, df2.iloc[i].values)[1])
                self.pValuesDf[f'Cluster {counter +1}'] = res
        elif test == 'kruskal':
            for counter, columns in enumerate(tqdm(self.sortedClusterList)):
                res = list()
                set2 = set(self.filterDf) - set(columns)
                df1, df2 = self.filterDf[columns], self.filterDf[set2]
                for i in tqdm(range(df1.shape[0])):
                    res.append(mannwhitneyu(df1.iloc[i].values, df2.iloc[i].values)[1])
                self.pValuesDf[f'Cluster {counter +1}'] = res
        #CORRECTED
        print('Performing FDR correction')
        for i in range(len(self.sortedClusterList)):
            self.pValuesDf['Cluster {} corr'.format(i+1)] = sm.stats.multitest.fdrcorrection(self.pValuesDf['Cluster {}'.format(i+1)].array, alpha=0.05, method='indep', is_sorted=False)[1]
        print('Filtering out genes')
        self.pValuesDf = self.pValuesDf.loc[(self.pValuesDf.iloc[:,4:] < 0.01).sum(axis = 1) < len(self.sortedClusterList)]
       # bla = pd.DataFrame(columns = [f'Cluster {i+1} corr:' for i in range(len(self.sortedClusterList))])
        # Sets the lesser significant genes to a p value of 1
       # tempDf = self.pValuesDf.iloc[:,len(self.sortedClusterList):].min(axis = 1)
       # for i in tqdm(range(self.pValuesDf.shape[0])):
       #     bla = bla.append(self.pValuesDf.iloc[:,len(self.sortedClusterList):].iloc[i].transform(lambda x: x if x == tempDf.iloc[i] else x/x))
       #self.pValuesDf.iloc[:,len(self.sortedClusterList):] = bla
        
            
    def filterEnrichedGenes(self, threshold = 2, filteringType = True, start = 1 ,transformed = True, csvpath = 'data/test.csv', split = '.', sep ='\t'):
        #Filters based on fpkm value. Not needed if using the test metho use threshold = 1
        if self.pValuesDf is not None:
            if filteringType == 'fpkm': self.readCsv(csvpath, readType = 'fpkm',  sep = sep, fpkm = True)
            elif filteringType == 'log2fold': self.logFoldChange(csvpath, sep = sep)
            #
            genes = list()
            #Fethes the significant genes
            for i in range(len(self.sortedClusterList)):
                genes.append(self.pValuesDf[self.pValuesDf['Cluster {} corr'.format(i+1)] < 0.01].index.tolist())
            self.tList = list()
            #Filters the genes based of fpkm
            if filteringType == 'fpkm':
                for i in range(len(self.sortedClusterList)):
                    finalList = genes[i]
                    lista = list()
                    for gene in finalList:
                        if  (self.fpkmDf[gene] > threshold):
                            lista.append(gene)
                    self.tList.append(lista)
            if filteringType == 'log2fold':
                for i in range(len(self.sortedClusterList)):
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
                    lista.append(gene.split(split,1)[0])
                self.eList.append(lista)
        else: print('Error: Run clusterTTest first')
            
    def log2filtering(self):
        #Filters groups with log2fold change lesser than 2
        self.gexDf = self.gexDf.filter(items = list(self.lfcDf[abs(self.lfcDf > 2)].dropna(axis=0,how='all').index), axis = 0)
        self.transposedGexDf = self.gexDf.transpose()

            
    def plotEnrichedGenes(self, cluster = 'all', clustered = True, vmin = -1, vmax = 1, row_cluster = False, **kwargs):
        #plots a clustermap of the working df
        if clustered: 
            if cluster == 'all':
                #Creates the dataframe to be plotted by filtering on the enriched gene set and stacking them
                #on top of each other cluster wise
                self.filterDf.index = [x.split('.')[0] for x in self.filterDf.index.tolist()]
                self.mapDf = pd.DataFrame(columns = list(self.gexDf))
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
                #     self.mapDf = pd.concat([self.mapDf, self.transposedGexDf.filter(items=genes)], axis = 1)
                ax = sns.clustermap(data = self.mapDf, cmap = 'vlag', metric = 'correlation', method = "average",vmax= vmax, vmin = vmin,col_cluster=False, col_colors=self.col_colors, row_colors = self.row_colors, row_cluster = row_cluster, **kwargs)
            else:
                sns.clustermap(data = self.filterDf.filter(items=self.tList[(cluster)], axis = 'index'), cmap = 'vlag', metric = 'correlation', method = "average",vmax= vmax, vmin = vmin,col_cluster=False, col_colors=self.col_colors, row_cluster = row_cluster, **kwargs)
        else: sns.clustermap(data = self.filterDf, cmap = 'vlag', metric = 'correlation', method = "average",vmax= vmax, vmin = vmin,col_cluster=False, col_colors=self.col_colors, row_cluster = row_cluster, **kwargs)
        
    
    def enrichR(self, secondRun = False, dbs = ['Enrichr_Libraries_Most_Popular_Genes']):
        #Fetches enrichment info from enrichR against different databases.
        #Fetches user list id for api calls
        if self.eList is not None:
            c = [0,1,2,3,4,5,6,7,8,'db', 'cluster']
            self.goEnrichDf = pd.DataFrame(columns = c)
            for counter, cluster in enumerate(self.eList):
                if self.samResults is not None:
                    if counter %2 == 0:
                        cNum = str(int(counter/2) + 1)+ ' up'
                    else:
                        cNum = str(int(counter/2) + 1)+ ' down'
                else:
                    cNum = counter
                print(f'Running for cluster {cNum}')
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
                    testDf['cluster'] = [cNum]*len(testDf.index.tolist())
                    self.goEnrichDf = self.goEnrichDf.append(testDf)
            self.goEnrichDf = self.goEnrichDf.sort_values(by = [2], ascending = False)
        else: print('Error: Run the filtering module')
            
    def plotEnrichmentResults(self,u_d = 'up',  **kwargs):
        #plots the enirchment result from rDf
        if self.goEnrichDf is not None:
            lista = self.goEnrichDf.sort_values('cluster')[1].tolist()
            meltDf = pd.melt(self.goEnrichDf, id_vars = ['cluster', 1], value_vars = [2])
            meltDf = meltDf.drop(columns = ['variable'])
            pivotDf = meltDf.pivot(index = 1, columns = 'cluster', values = 'value').transpose()
            pivotDf.index = set(lista)
            pal = sns.color_palette(n_colors = 6)
            col_colors = [pal[int(i/2)] for i in range(self.clusterNum*2)]
            col_colors = pd.Series(data = col_colors, index = pivotDf.columns.tolist())
            pivotDf = pivotDf.fillna(0)
            upcol = [col for col in pivotDf.columns.tolist() if 'up' in col]
            downcol = [col for col in pivotDf.columns.tolist() if 'down' in col]
            upDf = pivotDf[upcol]
            downDf = pivotDf[downcol]
            f = lambda x: x.max()
            upDf = upDf[(upDf.T != 0).any()]
            downDf = downDf[(downDf.T != 0).any()]
            downDf = pd.concat([downDf[f'{i+1} down'].sort_values(ascending=False)[downDf[f'{i+1} down'] > 0] for i in range(self.clusterNum)],axis = 'columns')
            upDf = pd.concat([upDf[f'{i+1} up'].sort_values(ascending=False)[upDf[f'{i+1} up'] > 0] for i in range(self.clusterNum)],axis = 'columns')
            
            downDf = downDf.fillna(0)
            upDf = upDf.fillna(0)
            
            if u_d == 'up':
                c  = sns.dark_palette("red")
            else:
                c  = sns.dark_palette("lightblue")
            dicts = {'up':upDf, 'down':downDf}
            #fig, ax = plt.subplots(figsize=figSize) 

            return sns.clustermap(dicts[u_d], vmax = 5, col_cluster = False,row_cluster = False, col_colors = col_colors, cmap = c,  **kwargs)
        else: print ('Run the enrichR module')
        
    def demoRun(self, cutoffs, csvpath, clusterpath, columnname = 'Gene', **kwargs):
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
                self.calcPcaModel(demo = cpath[-8:])
                self.Umap(3,demo = cpath[-8:])
                
        #plots the information from the fetched info
        g = sns.catplot(x= 'num', y='per_var',kind='bar', data = self.sDf, col='cutoff',col_wrap= 2, sharey = False, sharex = False, **kwargs)
        g2 = sns.relplot(x= 'comp 1', y='comp 2', data= self.uDf, style = 'hue', hue='hue', col='cutoff',col_wrap=2, sharey = False, sharex = False, **kwargs)
    
        
    def ensemble2gene(self, scopes = 'ensembl.gene'):
        #Converts the enrichment list from ensembl to gene symbol
        for i, ele in enumerate(self.eList):    
            self.eList[i] = self.ens2symHelper(ele, scopes )['symbol'].dropna().tolist()

    def metaGenes(self, **kwargs):
        if self.symbolDf is None: self.ensembl2symbol()
        if self.metaGeneDf is None: self.readCsv(readType = 'meta-gene')
        df = self.metaGeneDf
        df = df.iloc[:, 0:2]
        self.metaGenesBoxDf = pd.DataFrame(columns = ['Expression', 'Cluster', 'Metagene'])
        for x in df['Network cluster'].unique():
            for i, cluster in enumerate(self.sortedClusterList):
                values = (self.symbolDf[cluster].filter(items = list(df[df['Network cluster'] == x]['Gene Symbol']), axis = 0).mean(axis = 1).values)
                tempDf = pd.DataFrame({'Expression': values, 'Cluster':[i+1 for c in range(len(values))], 'Metagene':[x for c in range(len(values))]})
                self.metaGenesBoxDf = pd.concat([self.metaGenesBoxDf, tempDf])
        ax = sns.catplot(kind = 'box', data = self.metaGenesBoxDf, x = 'Cluster', y = 'Expression', col = 'Metagene', **kwargs)
    
        
        
    def ens2symHelper(self, genes, scopes = 'ensembl.gene'): 
        #Need to remove dots before before submitting gene list.
        headers = {'content-type': 'application/x-www-form-urlencoded'}
        params = f'q={",".join(genes)}&scopes={scopes}&fields=symbol'
        res = requests.post('http://mygene.info/v3/query', data=params, headers=headers)
        data = json.loads(res.text)
        return pd.json_normalize(data).set_index('query')
        
        
    def ensembl2symbol(self, scopes = 'ensembl.gene'):
        print('This may produce duplicate indexes which will be removed, keeping the first occurence')
        newIndex = [x.split('.')[0] for x in list(self.gexDf.index)]
        self.gexDf.index = newIndex
        self.transposedGexDf = self.gexDf.transpose()
        self.genes = self.ens2symHelper(newIndex, scopes)
        newIndex = (self.genes['symbol'].dropna())
        self.symbolDf = self.gexDf.reindex(list(newIndex.index))
        self.symbolDf.index = list(newIndex)
        self.symbolDf = self.symbolDf[~self.symbolDf.index.duplicated(keep='first')]
        self.tSymbolDf = self.symbolDf.transpose()
        self.tSymbolDf['Clusters'] = self.transposedGexDf['Clusters']

        
                
        
    def calcCentroids(self, centroidpath = 'data/extraData/wilkerson.2012.LAD.predictor.centroids.csv', ensembl = True, method = 'pearson', **kwargs):
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
                if self.cDf[cType].corr(self.ctDf[sample], method = method) > temp:
                    temp = self.cDf[cType].corr(self.ctDf[sample], method = method)
                    cEle = cType
            cList.append(cEle)
        #convert the list of cancer types to colors, then joining it with the col_colors
        cDict = {'TRU':'Red', 'PP':'Blue' , 'PI':'Yellow'}
        centroidClusters = [inner for outer in [[counter+1 for i in x] for counter, x in enumerate([x  for x in self.sortedClusterList]) ] for inner in outer]
        
        c2List = [cDict[x] for x in cList]
        self.centroid = pd.DataFrame(c2List, columns = ['Centroid'])
        self.centroid.index = self.transposedGexDf.index
        self.centroid['Clusters'] = centroidClusters
        self.centroid['Cancers'] = cList
        self.col_colors = pd.concat([self.col_colors, self.centroid['Centroid']], axis = 1)
        ax = sns.countplot(x="Clusters", hue="Cancers", data=self.centroid, hue_order =['TRU', 'PP','PI'], **kwargs)
                
        
    def logFoldChange(self, csvpath, sep = '\t'):
        self.readCsv(csvpath = csvpath, readType = 'lfc', sep = sep)
        self.lfcDf = pd.DataFrame()
        for counter, columns in enumerate(self.sortedClusterList):
            self.lfcDf = pd.concat([self.lfcDf, self.rawDf.filter(items = columns).mean(axis = 1)],axis  =1)
        self.lfcDf.columns = [f'Cluster {cluster}' for cluster in range(1, 1+len(self.sortedClusterList) )]
        tempDf = self.lfcDf.copy()
        for i, column in enumerate(list(self.lfcDf)):
            self.lfcDf[f'Cluster {i+1}'] = np.log2(tempDf.loc[:, tempDf.columns != f'Cluster {i+1}'].mean(axis=1).div(tempDf[f'Cluster {i+1}']).replace(0,1)).fillna(0)
            self.lfcDf[f'Cluster {i+1}'] = np.log2(tempDf[f'Cluster {i+1}'].div(tempDf.loc[:, tempDf.columns != f'Cluster {i+1}'].mean(axis=1)).replace(0,1)).fillna(0)
       
        self.lfcDf.index = [x.split('.')[0] for x in self.rawDf[self.rawDf.columns.tolist()[0]].tolist()]
        
        
    def normalityTest(self):
        pvals = list()
        for column in tqdm(list(self.transposedGexDf)):
            pvals.append(normaltest(self.transposedGexDf[column])[1])
        sns.displot(pvals, binwidth = 0.05)
        
    def kaplanMeier(self):
        if self.clinicalDf is None:
             self.readCsv(readType = 'Clinical')
        df = self.clinicalDf
        km = kmf() 
        
        for counter, cluster in enumerate(self.sortedClusterList):
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
                    try:
                        float(followup[i])
                        Died.append(1)
                        time.append(followup[i])
                    except:
                        continue
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
            self.blabla = km.event_table
            
    def coxHazard(self):
        if self.clinicalDf is None:
             self.readCsv(readType = 'Clinical')
        
        
    def addFeature(self, feature = None, attr = None, censor = None, title = None):
        if title is None: title = feature
        if self.clinicalDf is None:
             self.readCsv(readType = 'Clinical')
        patients = ['-'.join(x.split('-')[0:3]) for x in self.clusterList]
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
        
    def plotSingleGene(self, genes = list(), scopes = 'ensembl.gene', **kwargs):
        if self.symbolDf is None: self.ensembl2symbol(scopes)
        fig, axs = plt.subplots(len(genes), figsize = (10,10))
        for counter, gene in enumerate(genes):
            lista = [gene, 'Clusters']
            g  = sns.violinplot(data = self.tSymbolDf[lista], x='Clusters', y =gene, ax = axs[counter], **kwargs)
            g.axhline(0)
        fig.tight_layout()
            
    def plotMultipleGenes(self, geneList = list,scopes = 'ensembl.gene', **kwargs):
        if self.symbolDf is None: self.ensembl2symbol(scopes = 'ensembl.gene')
        if len(geneList) > 1: fig, axs = plt.subplots(len(geneList), figsize = (10,10))
        for counter, genes in enumerate(geneList):
            t = pd.melt(self.tSymbolDf.filter(items = genes+ ['Clusters']), id_vars = 'Clusters' ).drop(['variable'], axis = 1)
            if len(geneList) > 1: g =  sns.boxplot(data = t, x = 'Clusters', y = 'value', ax = axs[counter], **kwargs)
            else: g = sns.violinplot(data = t, x = 'Clusters', y = 'value', **kwargs)
            g.axhline(0)
        if len(geneList) > 1: fig.tight_layout()
            
    def readImmuneData(self):
        samples = ['-'.join(sample.split('-')[0:3]) for sample in self.gexDf.columns.tolist()]
        sortedL = [['-'.join(x.split('-')[0:3]) for x in cluster] for cluster in self.sortedClusterList]
        f = lambda l: ['-'.join(x.split('-')[0:3]) for x in l]
        fDot = lambda l: ['-'.join(x.split('.')[0:3]) for x in l]
        
        
        neoDf = pd.read_csv('data/immuneData/TCGA_PCA.mc3.v0.2.8.CONTROLLED.filtered.sample_neoantigens_10062017.tsv', sep = '\t')
        neoDf = neoDf[neoDf['sample'].isin(samples)]
        neoDf = neoDf.set_index('sample')
        
        pmhcDf = pd.read_csv('data/immuneData/TCGA_pMHC_SNV_sampleSummary_MC3_v0.2.8.CONTROLLED_170404.tsv', '\t')
        pmhcDf['barcode'] = f(pmhcDf['barcode'].tolist())
        pmhcDf = pmhcDf[pmhcDf['barcode'].isin(samples)]
        pmhcDf = pmhcDf.set_index('barcode')
        
        mutDf = pd.read_csv('data/immuneData/mutation-load_updated.txt', sep='\t')
        mutDf = mutDf[mutDf['Patient_ID'].isin(samples)]
        mutDf = mutDf.set_index('Patient_ID')
        
        
        mastDf = pd.read_csv('data/immuneData/TCGA_mastercalls.abs_tables_JSedit.fixed.txt',sep = '\t')
        mastDf['sample'] = f(mastDf['sample'].tolist())
        mastDf = mastDf[mastDf['sample'].isin(samples)]
        mastDf = mastDf.set_index('sample')
        
        absDf = pd.read_csv('data/immuneData/ABSOLUTE_scores.tsv', sep='\t')
        absDf['Unnamed: 0'] = f(absDf['Unnamed: 0'].tolist())
        absDf = absDf[absDf['Unnamed: 0'].isin(samples)]
        absDf = absDf.set_index('Unnamed: 0')
        
        leuDf = pd.read_csv('data/immuneData/TCGA_all_leuk_estimate.masked.20170107.tsv', sep ='\t')
        leuDf.columns = ['Type', 'Sample', 'Leukocyte Fraction']
        leuDf['Sample'] = f(leuDf['Sample'].tolist())
        leuDf = leuDf[leuDf['Sample'].isin(samples)]
        leuDf = leuDf.set_index('Sample')
        leuDf = leuDf[~leuDf.index.duplicated(keep='first')]
        
        cibDf = pd.read_csv('data/immuneData/TCGA.Kallisto.fullIDs.cibersort.relative.tsv', sep = '\t')
        cibDf['SampleID'] = fDot(cibDf['SampleID'].tolist())
        cibDf = cibDf[cibDf['SampleID'].isin(samples)]
        cibDf = cibDf[~cibDf['SampleID'].duplicated(keep='first')]
        cibDf = cibDf.set_index('SampleID')
        
        mitDf = pd.read_csv('data/immuneData/mitcr_sampleStatistics_20160714.tsv', sep='\t')
        mitDf = mitDf[mitDf['ParticipantBarcode'].isin(samples)]
        mitDf = mitDf.set_index('ParticipantBarcode')
        mitDf = mitDf[~mitDf.index.duplicated(keep='first')]
        
        subDf = pd.read_csv('data/immuneData/TCGASubtype.20170308.tsv', sep ='\t')
        subDf = subDf[subDf['pan.samplesID'].isin(samples)]
        subDf = subDf.set_index('pan.samplesID')
        
        hrdDf = pd.read_csv('data/immuneData/TCGA.HRD_withSampleID.txt',sep = '\t')
        hrdDf['sampleID'] = f(hrdDf['sampleID'].tolist())
        hrdDf = hrdDf[hrdDf['sampleID'].isin(samples)].set_index('sampleID')
        
        self.immuneDf = pd.concat([leuDf, cibDf, mastDf, absDf, pmhcDf, neoDf, hrdDf, subDf], axis = 1)
        temp = self.transposedGexDf['Clusters']
        temp.index = f(temp.index.tolist())
        self.immuneDf['Clusters'] = temp
        
    def corrDNA(self, **kwargs):
        if self.immuneDf is None: self.readImmuneData()
        f = lambda l: ['-'.join(s.split('-')[0:3]) for s in l]
        hDf = pd.DataFrame()
        for i, samples in enumerate(self.sortedClusterList):
            hDf[f'Cluster {i+1}'] = self.immuneDf.filter(items=f(samples), axis = 'index').corr()['Leukocyte Fraction']
        hDf = hDf.dropna()
        self.hDf = hDf.iloc[1:]
        g = sns.clustermap(self.hDf, cmap = 'vlag', yticklabels = True, row_cluster = False, **kwargs)

        
    def plotSignatures(self, **kwargs):
        df = pd.read_csv('data/extraData/signature_profile_sample_mSignatureDB.txt', sep = '\t')
        df['Tumor_Sample_Barcode'] = [x.replace('.', '-') for x in df['Tumor_Sample_Barcode'].tolist()]
        df = df[df['Tumor_Sample_Barcode'].isin(['-'.join(sample.split('-')[0:3]) for sample in self.samples])]
        newC = self.col_colors
        newC.index = ['-'.join(x.split('-')[0:3]) for x in list(newC.index)]
        self.sigDf = df.pivot(index='Signature', columns = 'Tumor_Sample_Barcode', values = 'Contribution')[self.samples]
        sns.clustermap(self.sigDf, col_colors = newC, row_cluster = True, col_cluster = False, **kwargs)
    
    def boxplotExternalData(self, cols = 6):
        if self.immuneDf is None: self.readImmuneData()
        df = self.immuneDf
        for column in list(df):
            if not is_numeric_dtype(df[column]): df = df.drop([column], axis = 1)
        df['id'] = df.index.tolist()
        self.newDf = pd.melt(df, id_vars = ['id', 'Clusters'])
        g = sns.catplot(data = self.newDf, kind = 'box', col = 'variable', y ='value', x = 'Clusters', col_wrap = cols, sharey = False)
        
        
    def samAnalysis(self, properties = None, spiral = False, dist = None, q = 5, fc = 2, expressionData = None):
        if input('Have you moved cluster folder to the expressionData folder (y/n) ') == 'y'.lower():
            cmd = f'cd SAMDEG/dist && java -jar SAMDEG.jar {properties} {spiral} {dist} {int(self.clusterNum)} {q} {fc} "{expressionData}"'
            ans = input(f'Following command will be run in os.system: {cmd} \n Do you want to continue? (y/n) ')
            if ans.lower() == 'y':
                os.system(cmd)
            else:
                print('Cancelling')
        else:
            print('Please do so first')
        
    def plotSamResults(self, vmin = -1, vmax = 1.5, resultsPath = '/Users/jacobkarlstrom/projekt/SRIQ/notebook/data/expressionData/SRIQ_10000itr_1200var_10r/10000/QC_Spiral(false)/Results_log_o_5/SRIQ_Data_in_5_ClusterOrder_ABS_Unique.txt',q = 4, lfc = 2, **kwargs):
        lfcPath = '/'.join(resultsPath.split('/')[0:-1])+'/SMSAM/1_Cluster_vs_All_SAM_T-test.txt'
        qlfcDf = pd.read_csv(lfcPath, sep = '\t')
        genes = qlfcDf[qlfcDf['q-value'] < q][abs(qlfcDf[qlfcDf['q-value'] < q]['FC']) > lfc].index.tolist()
        df = pd.read_csv(resultsPath, sep = '\t')
        cols = df.columns.tolist()[1:]
        
        df = df.iloc[:,:-1]
        df.columns = cols
        index = df.index.tolist()
        index = [i for i in index if i in genes]
        df = df.reindex(index)
        f = lambda x: x-x.median()
        #df = df.filter(items = genes, axis = 'index')
        df = df.transform(f, axis = 'columns')
        cols = self.gexDf.columns.tolist()
        df = df[cols]
        sns.clustermap(df, vmin = vmin, vmax  = vmax, cmap ='vlag', row_cluster = False, col_cluster = False, **kwargs)
        self.samResults = df
       
        resultsPath = '/'.join(resultsPath.split('/')[:-1]) + '/SMSAM/Genelist_ABS_q(5.0)_fc(2.0)'
        
        files = sorted(os.listdir(resultsPath))
        

        self.fcDf = pd.DataFrame(columns = ['Gene', 'FC', 'Clusters'])
        
        
        for counter, file in enumerate(files):
            df = pd.read_csv(resultsPath +'/'+ file, sep = '\t', header = None, names = ['Gene', 'FC'])
            df = df.sort_values('FC')
            df = pd.concat([df, pd.Series(name = 'Clusters', data = [counter +1]*df.shape[0])], axis = 1)
            self.fcDf = pd.concat([self.fcDf, df])
        
        
        self.tList = list()
        for i in range(1,self.clusterNum+1):
            self.tList.append(self.fcDf[self.fcDf['Clusters']==i][self.fcDf[self.fcDf['Clusters']==i]['FC'] < 0]['Gene'].tolist())
            self.tList.append(self.fcDf[self.fcDf['Clusters']==i][self.fcDf[self.fcDf['Clusters']==i]['FC'] > 0]['Gene'].tolist())
        self.eList = [[x.split('.')[0] for x in l] for l in self.tList]
        
        
        
        
        
        