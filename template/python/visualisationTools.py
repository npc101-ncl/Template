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

def GFID(myDict):
    """Get First in Dictionary
    
    A function to strip the first and assumed only eliment out of a
    dictionary.
    
    Args:
       myDict (Dict):  Dictioary wraped around single member.
       
    Returns:
       what ever was in the dictionary.
    """
    return (myDict[next(iter(myDict))])

def trim_axs(axs, N):
    """removes unwanted subplots from subplots function
    
    Convenience function used to remove unnneeded axises generated
    by the subplots function and return them in a flat structure
    
    Args:
       axs (array of Axes objects):  output axises from subplots
       N (int): number of subplots you actualy want.
       
    Returns:
       array of Axes objects: your remaining axises
    """
    axs = axs.flat
    for ax in axs[N:]:
        ax.remove()
    return axs[:N]

def getTCSelectionMaxes(TCList,selectionList=None,varSelection=None,
                       valRemove=[]):
    """Get TimeCourse maximum values by variable
    
    Gets the maximum value for each variable over all TimeCouses and
    all their time points. Useful in ordering and grouping multi plots.
    
    Args:
       TCList (list):  A list of TimeCourse Dataframes
       
    Kwargs:
       selectionList (list of ints): an optional list of the TimeCourses
           (by their indexing in the TCList) to be included in the
           calculation.
       varSelection (list of str): an optional list of variables from the
           TimeCourses we want returned in the results.
       valRemove (list of float): an optional list of values that are to be
           excluded from the list. A quick and dirty way to exclude
           variables that have all been set to some pre determined value.
       
    Returns:
       serise: the maximum value found for each variable
    """
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
    """breaks Series into subparts of similar size
    
    Takes a serise and breaks it into parts of comparable value size.
    Sub parts are returned sorted internaly and relative to each other in
    decending order. Useful for ordering and grouping plots with very
    difrent value ranges.
    
    Args:
       mySerise (serise):  A serise of floats to break into parts
       
    Kwargs:
       relativeScaleBreakPoint (float): Specifies the smallest (by scale)
           drop in value that should force an entry into a new subpart.
       maxRunLength (int): Max size of any subpart.
       
    Returns:
       list of list of str: the indexes of the serise are returned broken
       into subparts (the inner list) acording the the size of their
       asociated values.
    """
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
    
    def multiPlot2(self,indexSelect=None,varSelect=None,wrapNumber=5,
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
            
    def multiPlot(self,indexSelect=None,varSelect=None,wrapNumber=5,
                  compLines=None, save = None, xlim = None,
                  forceYAxisZero = True):
        """Plots grid if time course variables
        
        Creats grid of rainbow coloured plots for each variable in
        colection of TimeCourses. Optionaly plots refrence data against
        it as black dots.
        
        Kwargs:
           indexSelect (list of ints): Indexs of time courses that should
               be included in the plot.
           varSelect (list of str): Variables from the Time course that
               should be graphed.
           wrapNumber (int): max number of columns to alow in grid.
           compLines (dataFrame): refrence data to graph against time
               courses. Sould be in wide format and have the times as the
               index.
           save (str): A path to where to save the image. if omited image
               not saved
        """
        if compLines is not None:
            compVars=list(compLines.columns)
            dfB = compLines.copy()
            if "Time" not in compVars:
                dfB["Time"]=dfB.index
            else:
                compVars.remove("Time")
            dfB = pd.melt(dfB,id_vars=["Time"],
                          value_vars=compVars)
        if varSelect is None:
            varSelect=list(self.longData['variable'].unique())
        if indexSelect is None:
            indexSelect=list(self.longData['index'].unique())
        if not isinstance(indexSelect,list):
            indexSelect = [indexSelect]
        if len(varSelect)<wrapNumber:
            cols = math.floor(math.sqrt(len(varSelect)))
        else:
            cols = wrapNumber
        rows = math.ceil(len(varSelect)/cols)
        fig, axs = plt.subplots(rows, cols, sharex=True, figsize=(12,10))
        if (rows>1):
            axs = trim_axs(axs, len(varSelect))
        elif (cols==1):
            axs = [axs]
        myColorMap = plt.get_cmap(name="hsv", lut=len(indexSelect)+1)
        for ax, theVar in zip(axs, varSelect):
            ax.set_title(theVar)
            df = self.longData
            df = df[df['variable']==theVar]
            if indexSelect is not None:
                for theIndex, i in zip(indexSelect,range(len(indexSelect))):
                    df2 = df[df['index']==theIndex]
                    ax.plot(df2["Time"], df2["value"],linestyle='solid',
                            color=myColorMap(i))
            if compLines is not None:
                dfB2 = dfB[dfB['variable']==theVar]
                ax.plot(dfB2["Time"], dfB2["value"],"ko")
            if xlim is not None:
                ax.set_xlim(xlim)
            if forceYAxisZero:
                ax.set_ylim([0, None])
        fig.tight_layout()
        if save is not None:
            fig.savefig(save)
        
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
                   save = None, indexNames=None):
        """Creates Violin plot showing variation in paramiter space
        
        Kwargs:
           indexSelect (list of ints): Indexs of sets of paramiter
               estimation results that should be ploted. Takes a maximum of
               2 to allow comparision of two sets. optionaly a single
               interger may be used. Omision leads to all sets of paramiter
               estimations being processed in a combined way.
           varSelect (list of str): Variables from the paramiter
               estimations that should be graphed.
           RSSSelect (list of floats): Specifies a maximum RSS value for
               data in the paramiter estimations to be included in the plot.
               Ordering should match indexSelect.
           save (str): A path to where to save the image. if omited image
               not saved
           indexNames (list of str): names to use in the plot legend. If
               Omited numerical index is used. Ordering should match
               indexSelect.
        """
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
        if indexNames is not None:
            if indexSelect is not None:
                if len(indexSelect)!=len(indexNames):
                    return None
            df = df.copy()
            df["index"] = df["index"].replace(
                    dict(zip(indexSelect,indexNames)))
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
            
    def waterFall(self,save = None, indexNames=None):
        """creates watrfall plot of RSS values in paramiter estimations
        
        Kwargs:
           save (str): A path to where to save the image. if omited image
               not saved
           indexNames (list of str): names to use in the plot legend.
        """
        if indexNames is not None:
            if not isinstance(indexNames,list):
                indexNames = [indexNames]
            df = self.RSSData.copy()
            df["index"] = [indexNames[x] for x in df["index"]]
        else:
            df = self.RSSData
        plt.figure()
        lp = sns.lineplot(data=df,x="subIndex",y="RSS",hue="index")
        if save is not None:
            lp.get_figure().savefig(save)
            
    def refPointVsRSS(self, myRef, save=None, indexNames=None):
        names=list(myRef.index)
        refPoint = myRef.copy()
        refPoint["index"] = 1
        refPoint["RSS"] = 1
        diffTable = self.wideData.copy()
        diffTable = (diffTable/refPoint)-1
        diffTable["index"] = diffTable["index"]+1
        diffTable["index"] = diffTable["index"].astype(int)
        diffTable["RSS"] = diffTable["RSS"]+1
        diffTable["difrence"] = 0
        for name in names:
            diffTable["difrence"] = (diffTable["difrence"] +
                     diffTable[name]*diffTable[name])
        diffTable["difrence"] = diffTable["difrence"] / len(names)
        if indexNames is not None:
            if not isinstance(indexNames,list):
                indexNames = [indexNames]   
            diffTable["index"] = [indexNames[x] for x in diffTable["index"]]
        plt.figure()
        sp = sns.scatterplot(data=diffTable, x="RSS", y="difrence",
                             hue="index")
        if save is not None:
            sp.get_figure().savefig(save)