#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 16:39:30 2019

@author: peter
"""

from sklearn.cluster import DBSCAN, MeanShift, estimate_bandwidth
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

def RSSClusterEstimation(PEDataFrame):
    myScaledData = pp.StandardScaler().fit_transform(PEDataFrame.filter(
            items=["RSS"],))
    bandwidth = estimate_bandwidth(myScaledData)
    ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
    ms.fit(myScaledData)
    returnList = []
    for i in range(len(set(ms.labels_))):
        returnList.append({"id":i, "maxRSS":0, "size":0})
    for i in range(len(ms.labels_)):
        returnList[ms.labels_[i]]["maxRSS"] = max(
                returnList[ms.labels_[i]]["maxRSS"],
                PEDataFrame["RSS"].iloc[i])
        returnList[ms.labels_[i]]["size"]=returnList[ms.labels_[i]]["size"]+1
    returnList.sort(key=lambda cl:cl["maxRSS"])
    return returnList