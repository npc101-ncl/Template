#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 12:03:49 2021

@author: peter
"""

from python.visualisationTools import *
from python.analysisTools import *
from python.utilityTools import *
import os, re, sys
import pandas as pd
import pickle
import seaborn as sns
import matplotlib.pyplot as plt

name = "Hill_coop"

working_directory = os.path.dirname(os.path.abspath(__file__))
RS = loadPick(["data",name,'runSwitches.p'],relative=True)
cols=5

for myType, mySet in RS["PLSets"].items():
    for paramCase in mySet:
        myPL = loadPick(["data",name,"proLick-"+str(paramCase)+".p"],
                     relative=True)
        df = []
        for variable in myPL[0].keys():
            tdf = myPL[0][variable][[variable,"RSS"]].copy()
            tdf= tdf.rename(columns={variable:"adjustment"})
            tdf["variable"] = variable
            if variable in myPL[1].keys():
                if myPL[1][variable]!=0:
                    tdf["adjustment"] = tdf["adjustment"]/myPL[1][variable]
                else:
                    print(variable,":",myPL[1][variable])
            else:
                print("missing",variable)
            df.append(tdf)
        df = pd.concat(df, ignore_index=True)
        rows = len(myPL[0].keys())//cols
        if len(myPL[0].keys())%cols>0:
            rows=rows+1
        fig, axs = plt.subplots(rows, cols, sharex=True, figsize=(12,10))
        axs = trim_axs(axs,len(myPL[0].keys()))
        for ax, variable in zip(axs,myPL[0].keys()):
            ax.set_xscale('log')
            ax.title.set_text(variable)
            ax.plot(df[df["variable"]==variable]["adjustment"],
                    df[df["variable"]==variable]["RSS"])
        fig.tight_layout()
        os.makedirs(resolvePath(['figures', name],relative=True),exist_ok=True)
        fig.savefig(os.path.join(working_directory,'figures', name, "PL"+
                                 str(paramCase)+"-"+myType+".png"))
