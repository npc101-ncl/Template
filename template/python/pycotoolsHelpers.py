#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 14:00:40 2019

@author: peter
"""

import site, os, time, re
from pycotools3 import model, tasks, viz
import pandas as pd

def addCopasiPath(copasiPath):
    if not re.search(copasiPath, os.environ["PATH"]):
        os.environ["PATH"] += os.pathsep + copasiPath
        
def extractAntReactions(antString):
    antLines=antString.splitlines()
    reactions = []
    for antLine in antLines:
        reaction={}
        temp=antLine.split("->")
        if len(temp) != 2:
            temp=antLine.split("=>")
            if len(temp) != 2:
                continue
            else:
                reaction['irreversible']=True
        else:
            reaction['irreversible']=False
        temp2=temp[0].split(":")
        if len(temp2)!=2:
            if len(temp2)==1:
                LHS=temp2[0].strip()
            else:
                continue
        else:
            reaction['name']=temp2[0].strip()
            LHS=temp2[1].strip()
        temp2=temp[1].split(";")
        if len(temp2)==1:
            RHS=temp2[0].strip()
        else:
            reaction['formula']=temp2[1].strip()
            RHS=temp2[0].strip()
        LHS=[i.strip() for i in LHS.split("+")]
        if LHS==['']:
            LHS=[]
        for i in range(len(LHS)):
            if LHS[i][0]=='$':
                LHS[i]={'fixed':True, 'var':LHS[i][1:]}
            else:
                LHS[i]={'fixed':False, 'var':LHS[i]}
        reaction['LHS']=LHS
        RHS=[i.strip() for i in RHS.split("+")]
        if RHS==['']:
            RHS=[]
        for i in range(len(RHS)):
            if RHS[i][0]=='$':
                RHS[i]={'fixed':True, 'var':RHS[i][1:]}
            else:
                RHS[i]={'fixed':False, 'var':RHS[i]}
        reaction['RHS']=RHS
        reactions.append(reaction)
    return reactions

def replaceVariable(theString,oldName,newName):
    return newName.join(re.split('(?<![\w_])'+oldName+'(?![\w_])',
                                 theString))
    
def renameCSVColumns(CSVPathOrList,targetDir,renameDict):
    if not isinstance(CSVPathOrList, list):
        CSVPathOrList = [CSVPathOrList]
    outList=[]
    for CSVPath in CSVPathOrList:
        CSVData = pd.read_csv(CSVPath)
        CSVData=CSVData.rename(columns=renameDict)
        newPath = os.path.join(targetDir,os.path.split(CSVPath)[1])
        CSVData.to_csv(newPath,index=False)
        outList.append(newPath)
    if len(outList)==1:
        outList=outList[0]
    return outList
            
  
class modelRunner:
    def __init__(self, antString=None, run_dir=None):
        if run_dir is None:
            # run_dir = os.path.dirname(os.path.abspath(''))
            run_dir = os.path.abspath('')
            run_dir = os.path.join(run_dir, 'copasiRuns')
            i=1
            nameFree=False
            while not nameFree:
                temp_dir = os.path.join(run_dir, 'run'+str(i))
                if not os.path.isdir(temp_dir):
                    os.makedirs(temp_dir)
                    self.run_dir=temp_dir
                    nameFree=True
                else:
                    i=i+1
        else:
            if os.path.isdir(run_dir):
                self.run_dir=run_dir
            else:
                os.makedirs(run_dir)
                self.run_dir=run_dir
        if antString is None:
            self.antString = '''
            model empty_model()
            end
            '''
        else:
            self.antString=antString
            
    def genPrefixAntString(self,estimatedVar,prefix="_"):
        self.prefixAntString = self.antString
        for name in estimatedVar:
            self.prefixAntString = replaceVariable(self.prefixAntString,
                                                   name,prefix+name)

    def runParamiterEstimation(self,expDataFP,PEName=None,
                               copyNum=1,estimateIC=True,estimatedVar=None,
                               prefix='_',rocket=False):
        if PEName is None:
            copasi_filename = os.path.join(self.run_dir,
                                           'paramiterEstimation.cps')
        else:
            copasi_filename = os.path.join(self.run_dir, PEName)
        
        if estimatedVar is not None:
            self.genPrefixAntString(estimatedVar)
            self.recentModel = model.loada(self.prefixAntString,
                                           copasi_filename)
            colRenameDict = {i:(prefix+i) for i in estimatedVar}
            expDataFP = renameCSVColumns(expDataFP,self.run_dir,
                                         colRenameDict)
        else:
            self.recentModel = model.loada(self.antString, copasi_filename)
        
        if estimateIC:
            PEParams='gm'
        else:
            PEParams='g'
            
        if rocket==True:
            runMode='slurm'
        else:
            runMode=True
        
        with tasks.ParameterEstimation.Context(self.recentModel, expDataFP, 
                                               context='s', 
                                               parameters=PEParams) as context:
            context.set('randomize_start_values', True)
            context.set('separator', ',')
            if estimatedVar is not None:
                context.set('prefix', prefix)
            context.set('method', 'genetic_algorithm')
            context.set('population_size', 50)
            context.set('number_of_generations', 300)
            context.set('run_mode', runMode) 
            context.set('pe_number', 1) 
            context.set('copy_number', copyNum) 
            config = context.get_config()
            
        self.recentPE = tasks.ParameterEstimation(config)
        doneNum=0
        while copyNum>doneNum:
            try:
                parse_object = viz.Parse(self.recentPE)
                if isinstance(parse_object[list(parse_object)[0]],
                                          pd.DataFrame):
                    doneNum = parse_object[list(parse_object)[0]].shape[0]
            except:
                time.sleep(10)
            if rocket:
                print(str(doneNum) + ' of ' + str(copyNum) + ' done')
        
        if estimatedVar is not None:
            colRenameDict = {(prefix+i):i for i in estimatedVar}
        return_data={}
        for key in list(parse_object):
            if estimatedVar is None:
                return_data[key]=parse_object[key].copy()
            else:
                return_data[key]=parse_object[key].copy()
                return_data[key].rename(columns=colRenameDict, inplace=True)
        return return_data
    
    def runTimeCourse(self,duration,stepSize=0.01,intervals=100,
                      PEName=None,adjustParams=None,subSet=None):
        if PEName is None:
            copasi_filename = os.path.join(self.run_dir,
                                           'timeCourse.cps')
        else:
            copasi_filename = os.path.join(self.run_dir, PEName)
            
        self.recentModel = model.loada(self.antString, copasi_filename)
        if adjustParams is not None:
            if subSet is None:
                subSet = range(len(adjustParams.index))
            results = []
            for setIndex in subSet:
                model.InsertParameters(self.recentModel,df=adjustParams,
                                       index=setIndex,inplace=True)
                self.recentTimeCourse = tasks.TimeCourse(self.recentModel,
                                                         end=duration,
                                                         step_size=stepSize,
                                                         intervals=intervals)
                results.append(viz.Parse(self.recentTimeCourse).data.copy())
            return results
        else:
            self.recentTimeCourse = tasks.TimeCourse(self.recentModel,
                                                     end=duration,
                                                     step_size=stepSize,
                                                     intervals=intervals)
            return viz.Parse(self.recentTimeCourse).data.copy()