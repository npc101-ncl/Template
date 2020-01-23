#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 14 01:00:04 2019

@author: peter
"""

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import math

sns.set(context='paper')

def getTCSelectionMaxes(TCList,selectionList=None,varSelection=None,
                       valRemove=[]):
    if selectionList is None:
        maxes = [course.max() for course in TCList]
    else:
        maxes = [TCList[selection].max() for selection in selectionList]
    maxes = pd.DataFrame(maxes).max()
    maxes = maxes.drop(labels=['Time'])
    if varSelection is not None:
        maxes = maxes.filter(items=varSelection)
    maxes = maxes.sort_values(ascending=False)
    for value in valRemove:
        maxes = maxes[maxes!=value]
    return maxes

def breakSeriesByScale(mySerise,relativeScaleBreakPoint=0.1,
                       maxRunLength=6):
    outerList=[]
    innerList=[]
    lastVal=None
    for index, value in mySerise.sort_values(ascending=False).items():
        if lastVal is None:
            innerList.append(index)
            lastVal=value
            continue
        if (lastVal*relativeScaleBreakPoint>value or 
            len(innerList)==maxRunLength):
            outerList.append(innerList)
            innerList=[]
        else:
            innerList.append(index)
        lastVal=value
        outerList = [myList for myList in outerList if len(myList)>=1]
    return outerList

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
    
    def multiPlot(self,indexSelect=None,varSelect=None,wrapNumber=5,
                  compLines=None, save = None):
        df = self.longData
        if varSelect is not None:
            if not isinstance(varSelect, list):
                varSelect = [varSelect]
            df = df[df['variable'].isin(varSelect)]
        if indexSelect is not None:
            if not isinstance(indexSelect, list):
                indexSelect = [indexSelect]
            df=df.loc[df['index'].isin(indexSelect)]
        indexes = list(df['index'].unique())
        colors = sns.husl_palette(len(indexes)).as_hex()
        if compLines is not None:
            compVars=list(compLines.columns)
            dfB = compLines.copy()
            dfB["Time"]=dfB.index
            dfB = pd.melt(dfB,id_vars=["Time"],
                          value_vars=compVars)
            dfB["index"]=-1
            df = pd.concat([dfB,df],ignore_index=True)
            colors.append("#000000")
            indexes.append(-1)
        #colors = dict(zip(indexes, colors))
        grid = sns.FacetGrid(df, col="variable", col_wrap=wrapNumber, 
                             palette=colors)
        grid.map(sns.lineplot,"Time","value","index")
        if save is not None:
            grid.savefig(save)
        
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
            if "" in theColumns:
                df=df.drop(columns="")
            theColumns=list(df.columns).remove('RSS')
            df['subIndex'] = df.index
            df=pd.melt(df,id_vars=['subIndex','RSS'],
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
        
    def violinPlot(self,indexSelect=None,paramSelect=None,RSSSelect=None,
                   save = None):
        df=self.BPData
        if RSSSelect is not None:
            if not isinstance(RSSSelect, list):
                RSSSelect = [RSSSelect]
            if indexSelect is None:
                df = df.loc[df['RSS']<=RSSSelect[0]]
        if indexSelect is not None:
            if not isinstance(indexSelect, list):
                indexSelect = [indexSelect]
            if RSSSelect is not None:
                for i in range(len(indexSelect)):
                    df = df.loc[(df['RSS']<=RSSSelect[i]) |
                            (df['index']!=indexSelect[i])]
            df = df.loc[df['index'].isin(indexSelect)]
        if paramSelect is not None:
            if not isinstance(paramSelect, list):
                paramSelect = [paramSelect]
            df = df.loc[df['variable'].isin(paramSelect)]
        plt.figure()
        if indexSelect is None:
            if paramSelect is None:
                vp = sns.violinplot(x="variable", y="value", data=df, cut=0)
            else:
                vp = sns.violinplot(x="variable", y="value", data=df,
                                    order=paramSelect, cut=0)
        else:
            if paramSelect is None:
                vp = sns.violinplot(x="variable", y="value", data=df,
                                    hue="index", split=True, cut=0)
            else:
                vp = sns.violinplot(x="variable", y="value", data=df,
                                    hue="index", split=True,
                                    order=paramSelect, cut=0)
        if save is not None:
            vp.get_figure().savefig(save)
            
    def waterFall(self,save = None):
        plt.figure()
        lp = sns.lineplot(data=self.RSSData,x="subIndex",y="RSS",hue="index")
        if save is not None:
            lp.get_figure().savefig(save)