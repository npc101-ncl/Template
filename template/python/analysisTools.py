#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 16:39:30 2019

@author: peter
"""

from sklearn.cluster import DBSCAN
from sklearn import preprocessing as pp
import numpy as np
import pandas as pd

def clusterParameterEstimation(PEDataFrame,testDistance=5,minInCluster=2):
    myScaledData = pp.StandardScaler().fit_transform(PEDataFrame.drop(
            columns="RSS"))
    clustering = DBSCAN(eps=testDistance,
                        min_samples=minInCluster).fit(myScaledData)
    returnDict = {}
    for i in set(clustering.labels_):
        returnDict[str(i)]=[]
    for i in range(len(clustering.labels_)):
        returnDict[str(clustering.labels_[i])].append(i)
    return returnDict