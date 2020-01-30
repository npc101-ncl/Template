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
    """Clusters estimated paramiters
    
    based on their 'nearness' in paramiter space using DBSCAN

    Args:
       PEDataFrame (DataFrame):  paramiter estimations to cluster

    Kwargs:
       testDistance (float): eps value used in DBSCAN. Basicly the
           minimum 'distance' alowed between clusters.
       minInCluster (int): minimum size alowed for cluster.

    Returns:
       Dict of Lists: Each cluster is returned as an intiger indexed
       List except for the unclustered points which are returned as a
       list in the -1 dict entry.
    """
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
    """Clustered paramiter estimationes
    
    bassed purely on their RSS values using MeanShift.
    
    Args:
       PEDataFrame (DataFrame):  paramiter estimations to cluster
       
    Returns:
       Lists of Dicts: Clusteres are returned as dictionaries in order of
       the size of their largest RSS value (assending). dict structure is::
           
           {"id":arbitrary_id_assigned_to_cluster,
           "maxRSS":largest_RSS_value_in_cluster,
           "size":number_of_members_of_cluster}
           
    """
    myScaledData = pp.StandardScaler().fit_transform(PEDataFrame.filter(
            items=["RSS"]))
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