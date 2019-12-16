#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 14 01:00:04 2019

@author: peter
"""

import seaborn as sns
import pandas as pd

sns.set(context='paper')

class timeCourseVisualiser:
    def __init__(self,data):
        if not isinstance(data, list):
            data = [data]
        runningList=[]
        for i in range(len(data)):
            valueColumns=list(data[i].columns)
            if "Time" in valueColumns:
                valueColumns.remove("Time")
            if "" in valueColumns:
                valueColumns.remove("")
            dataFrame=pd.melt(data[i],id_vars=["Time"],
                    value_vars=valueColumns)
            dataFrame['index']=i
            runningList.append(dataFrame)
        self.longData = pd.concat(runningList,ignore_index=True)
    
    def multiPlot(self,indexSelect=None):
        grid = sns.FacetGrid(self.longData, col="variable")
        if indexSelect is None:
            grid.map(sns.lineplot,data=self.longData,x="Time",y="value")
        else:
            if not isinstance(indexSelect, list):
                indexSelect = [indexSelect]
            df=self.longData.loc[self.longData['index'].isin(indexSelect)]
            grid.map(sns.lineplot,data=self.longData,x="Time",y="value",
                     hue="index")
        
class parameterEstimationVisualiser:
    def __init__(self,data):
        if not isinstance(data, list):
            data = [data]
        BPList=[]
        RSSList=[]
        wideList=[]
        for i in range(len(data)):
            df = data[i][list(data[i])[0]].copy()
            theColumns=list(df.columns)
            if "RSS" in theColumns:
                df=df.drop(columns="RSS")
            if "" in theColumns:
                df=df.drop(columns="")
            theColumns=list(df.columns)
            df['subIndex'] = df.index
            df=pd.melt(df,id_vars=['subIndex'],
                       value_vars=theColumns)
            df['index']=i
            BPList.append(df)
            df = data[i][list(data[i])[0]].copy()
            df['subIndex'] = df.index
            df=df[['subIndex',"RSS"]]
            df['index']=i
            RSSList.append(df)
            df = data[i][list(data[i])[0]].copy()
            df['index']=i
            wideList.append(df)
        self.BPData = pd.concat(BPList,ignore_index=True)
        self.RSSData = pd.concat(RSSList,ignore_index=True)
        self.wideData = pd.concat(wideList,ignore_index=True)
        
    def violinPlot(self,indexSelect=None,paramSelect=None):
        df=self.BPData
        if indexSelect is not None:
            if not isinstance(indexSelect, list):
                indexSelect = [indexSelect]
            df = df.loc[df['index'].isin(indexSelect)]
        if paramSelect is not None:
            if not isinstance(paramSelect, list):
                paramSelect = [paramSelect]
            df = df.loc[df['variable'].isin(paramSelect)]
        if indexSelect is None:
            if paramSelect is None:
                sns.violinplot(x="variable", y="value", data=df, cut=0)
            else:
                sns.violinplot(x="variable", y="value", data=df,
                               order=paramSelect, cut=0)
        else:
            if paramSelect is None:
                sns.violinplot(x="variable", y="value", data=df,
                               hue="index", split=True, cut=0)
            else:
                sns.violinplot(x="variable", y="value", data=df,
                               hue="index", split=True, order=paramSelect,
                               cut=0)
            
    def waterFall(self):
        sns.lineplot(data=self.RSSData,x="subIndex",y="RSS",hue="index")