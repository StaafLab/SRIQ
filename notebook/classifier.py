#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 13:29:51 2020

@author: jacob
"""

import numpy as np
from sklearn.model_selection import train_test_split
from sklearn import datasets
from sklearn import svm
import numpy as np
import pandas as pd
from sklearn.model_selection import cross_val_score



class classifier():
    def __init__(self, df):
        self.df = df
    
    def divideDf(self, factorList):
        """" 
        Purpose: Divides a dataframe to an evenly sized dataframe
        Parameters: factorsList: list of factors that you want to base the division on
        Returns: X_train, X_validation, Y_train, Y_validation
        """
        factors1 = dict()
        nDf = self.df.transpose().filter(items=factorList).transpose()
        for column in list(nDf.columns):
            if tuple(nDf[column].values) in factors1:
                factors1[tuple(nDf[column].values)] = factors1[tuple(nDf[column].values)]+1
            else:
                factors1[tuple(nDf[column].values)] = 1
        
        retList1, retList2 = list(nDf.columns), list()
        factors2 = dict()
        flag = 0
        for counter, column in enumerate(list(nDf.columns)):
            if tuple(nDf[column].values) not in factors2:
                factors2[tuple(nDf[column].values)] = 1
                factors1[tuple(nDf[column].values)] = factors1[tuple(nDf[column].values)] -1
                retList2.append(retList1.pop(counter-flag))
                flag +=1
                
            elif factors1[tuple(nDf[column].values)] > factors2[tuple(nDf[column].values)]:
                factors1[tuple(nDf[column].values)] = factors1[tuple(nDf[column].values)] -1
                factors2[tuple(nDf[column].values)] = factors2[tuple(nDf[column].values)] +1
                retList2.append(retList1.pop(counter-flag))
                flag += 1
        return self.df.filter(items = retList1).transpose().drop(columns = ['Clusters']), self.df.filter(items = retList2).transpose().drop(columns=['Clusters']), self.df.filter(items = retList1).transpose()['Clusters'], self.df.filter(items = retList2).transpose()['Clusters']
        
        
        
        
    def SVC(self, validation='cross', factors = ['Clusters']):
        if validation == 'cross':
            X = self.df.transpose().drop(columns=['Clusters'])
            y = self.df.transpose()['Clusters']
            clf = svm.SVC(kernel='linear', C=1)
            scores = cross_val_score(clf, X, y, cv=10)
            return scores
        if validation == 'div':
            X_train, X_test, y_train, y_test = self.divideDf(factors)
            clf = svm.SVC(kernel='linear', C=1).fit(X_train, y_train)
            return clf.score(X_test, y_test)
            