#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 14:50:41 2021

@author: peter
"""

from python.utilityTools import *
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import math

def dfPreper(myTC,selectors,sigIden,gridAxis,values):
    temp = myTC
    for key, value in selectors.items():
        temp = [tc for tc in temp if tc.iloc[0][key]==value]
    for i in range(len(temp)):
        for j, myId in enumerate(sigIden):
            for key, value in myId.items():
                if abs(temp[i].iloc[0][key]-value)>0.0001:
                    break
            else:
                temp[i]["sigID"]=j
                break
    temp = pd.concat(temp, ignore_index=True)
    temp = pd.melt(temp, id_vars=["Time","sigID"]+gridAxis,
                   value_vars=values)
    temp2 = temp.groupby(by=gridAxis+["sigID","variable","Time"])
    temp3 = []
    temp3.append(temp2.mean())
    temp3[0] = temp3[0].rename(columns={"value": "mean"}, errors="raise")
    temp3.append(temp2.min())
    temp3[1] = temp3[1].rename(columns={"value": "min"}, errors="raise")
    temp3.append(temp2.max())
    temp3[2] = temp3[2].rename(columns={"value": "max"}, errors="raise")
    temp3 = pd.concat(temp3,axis=1)
    temp3.reset_index(inplace=True)
    return temp3

def getMinRow(df,minRow=None):
    mySerise = df[df[minRow]==df[minRow].min()].iloc[0]
    mySerise["row"] = mySerise.name
    return mySerise

cmdDict, cmdFlags = getCmdLineArgs()

if "name" in cmdDict.keys():
    name = cmdDict["name"]
else:
    name = "Hill_coop"

params = loadPick(["data",name,"eqParams.p"], relative=True)
RS = loadPick(["data",name,"runSwitches.p"], relative=True)
myTC = loadPick(["data",name,"timeCourses.p"], relative=True)

params = GFID(params)

for index, row in RS["iParamsDF"].iterrows():
    temp = params
    for col, item in row.items():
        temp = temp[temp[col]==item]
    print(temp[RS["oParams"][:5]].std())
    print(temp[RS["oParams"][5:]].std())
    
myStd = params.groupby(list(RS['iParamRanges'])).std()

print(myStd.mean()[RS["oParams"][:5]].mean(),
      myStd.mean()[RS["oParams"][5:]].mean())

resample = len(myTC)//len(RS['signalSet'])
TCG = [pd.concat(myTC[i:len(myTC):resample],ignore_index=True) for i
       in range(resample)]

for i in range(len(TCG)):
    TCG[i]["pro2v1"]=abs(TCG[i]["delay2A"]/TCG[i]["delay1"]-1)
    TCG[i]["pro3v1"]=abs(TCG[i]["delay3A"]/TCG[i]["delay1"]-1)
TCSum = pd.DataFrame([TC.mean() for TC in TCG])

badCandidates = []
goodCandidates = []
medianCandidates = []

myGroup = TCSum.groupby(list(RS['iParamRanges'].keys()))
temp = myGroup.apply(getMinRow, minRow="pro2v1")
temp = temp.reset_index(drop=True)
temp = temp.sort_values(by=['pro2v1'],ascending=False)
temp[list(RS['iParamRanges'].keys())+['pro2v1']]
for k in RS['iParamRanges'].keys():
    fig = plt.figure()
    s = plt.scatter(temp[k], temp['pro2v1'])
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel(k)
    plt.ylabel('pro2v1')
    plt.ylim(min(temp['pro2v1']), )
    
badCandidates.append(int(temp.iloc[0]["row"]))
medianCandidates.append(int(temp.iloc[len(temp)//2]["row"]))
goodCandidates.append(int(temp.iloc[-1]["row"]))


temp = myGroup.apply(getMinRow, minRow="pro3v1")
temp = temp.reset_index(drop=True)
temp = temp.sort_values(by=['pro3v1'],ascending=False)
temp[list(RS['iParamRanges'].keys())+['pro3v1']]
for k in RS['iParamRanges'].keys():
    fig = plt.figure()
    s = plt.scatter(temp[k], temp['pro3v1'])
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel(k)
    plt.ylabel('pro3v1')
    plt.ylim(min(temp['pro3v1']), )
badCandidates.append(int(temp.iloc[0]["row"]))
medianCandidates.append(int(temp.iloc[len(temp)//2]["row"]))
goodCandidates.append(int(temp.iloc[-1]["row"]))

def myPlot(candidates):
    for row in candidates:
        fig = plt.figure()
        for i, colour in zip(range(len(RS['signalSet'])),
                             ['b','g','r','c','m','y','k']):
            span = len(TCG[row])//len(RS['signalSet'])
            chunk = TCG[row].iloc[i*span:(i+1)*span]
            for style, output in zip(["-",":","--"],
                                      ['delay1','delay2A','delay3A']):
                plt.plot(chunk["Time"], chunk[output],colour+style)
        myT = TCSum.iloc[row][list(RS['iParamRanges'].keys())].to_dict()
        myT = {k:round(v,3) for k,v in myT.items()}
        plt.title(str(myT))
        
myPlot(badCandidates)
myPlot(medianCandidates)
myPlot(goodCandidates)

RS["PLSets"]={"badCandidates":badCandidates,
              "medianCandidates":medianCandidates,
              "goodCandidates":goodCandidates}
savePick(["data",name,"runSwitches.p"], RS, relative=True)

"""
gridAxis = ["kD1_S","kD1_V"]
temp = {k:len(v) for k,v in RS['iParamRanges'].items()
        if not k in gridAxis}
myLen = 1
for _,i in temp.items():
    myLen *= i
for i in range(myLen):
    selectors={}
    j=i
    for k,v in temp.items():
        selectors[k]=RS['iParamRanges'][k][j%v]
        j=j//v
    fig, axs = plt.subplots(len(RS['iParamRanges'][gridAxis[0]]), 
                            len(RS['iParamRanges'][gridAxis[1]]),
                            sharex=True, sharey=True)
    for i, row in enumerate(axs):
        for j, cell in enumerate(row):
            subSelectors = selectors
            subSelectors[gridAxis[0]] = RS['iParamRanges'][gridAxis[0]][i]
            subSelectors[gridAxis[1]] = RS['iParamRanges'][gridAxis[1]][j]
            df = dfPreper(myTC, subSelectors, RS['signalSet'], [],
                          ['delay1','delay2A','delay3A'])
            for style, variable in zip(["-",":","--"],
                                       ['delay1','delay2A','delay3A']):
                for hue, sigID in zip(['b','g','r','c','m','y','k'],
                                      range(len(RS['signalSet']))):
                    df2=df[df["sigID"]==sigID]
                    df2=df2[df2["variable"]==variable]
                    cell.fill_between(df2["Time"],df2["min"],df2["max"],
                                      ls=style,color=hue,alpha=0.5)
"""