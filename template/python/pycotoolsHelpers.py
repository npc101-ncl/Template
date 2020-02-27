#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 14:00:40 2019

@author: peter
"""

import site, os, time, re
from pycotools3 import model, tasks, viz
import tellurium as te
import pandas as pd
import logging

def addCopasiPath(copasiPath):
    """adds the path to copasi to pythons working path
    
    Copasi is often not in the path python uses. This will check if its
    not and add it.

    Args:
       copasiPath (str):  a string with the path to copasi on the
       relivent machine.
    """
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
    """changes the name of a variable in an antimony string
    
    It's helpfull to be able to rename variables in antimony strings on
    the fly to utalise pycotools 'prefix' feature. This quick a dirty
    methiod uses regular expresions to do this and works suprisingly
    reliably.

    Args:
       theString (str):  the antimony string to parse
       oldName (str):  name of variable to replace
       newName (str):  name of variable to be the replacement

    Returns:
       str: the new modified antimony string.
    """
    return newName.join(re.split('(?<![\w_])'+oldName+'(?![\w_])',
                                 theString))
    
def renameCSVColumns(CSVPathOrList,targetDir,renameDict,indepToAdd=None):
    """preproceses CSV data files for paramiter estimations
    
    Allows us to take one of more CSV files and rename their columns (for
    example to take advndage of prefix feture) and optionaly add new
    columns bearing the "_indep" sufix used to indicate a column should
    be treated as an indipendent variable in the paramiter estimation.
    
    Args:
       CSVPathOrList (str or list of str):  path or list of paths to csv
           files to modify.
       targetDir (str):  path to directory in which modified csv files
           should be saved.
       renameDict (dict):  dictioary of str key value pairs that indicate
           the replacments to perform ineach csv files column names.
       
    Kwargs:
       indepToAdd (dict or list of dicts): key value pairs representing
           the independed variables that should be added to each data set
           and the corisponding values they should be set to. Ordering of 
           list should match CSVPathOrList.
           
    Returns:
       List or str: path or list of paths to saved modified files.
       (ordering will match CSVPathOrList)
    """
    if not isinstance(CSVPathOrList, list):
        CSVPathOrList = [CSVPathOrList]
    if indepToAdd is not None:
        if not isinstance(indepToAdd, list):
            indepToAdd = [indepToAdd]
        if len(CSVPathOrList)!=len(indepToAdd):
            print("each data set needs independed variables if you use them")
            return None
    outList=[]
    for CSVPath, i in zip(CSVPathOrList,range(len(CSVPathOrList))):
        CSVData = pd.read_csv(CSVPath)
        CSVData=CSVData.rename(columns=renameDict)
        if indepToAdd is not None:
            for myKey, myVal in indepToAdd[i].items():
                CSVData[myKey+"_indep"]=myVal
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
        self.methodDict = {
                "particle_swarm_default":{
                        "method":"particle_swarm",
                        "swarm_size":50,
                        "iteration_limit":2000},
                "particle_swarm_rigorous":{
                        "method":"particle_swarm",
                        "swarm_size":100,
                        "iteration_limit":4000},
                "particle_swarm_aggressive":{
                        "method":"particle_swarm",
                        "swarm_size":200,
                        "iteration_limit":6000},
                "particle_swarm_heroic":{
                        "method":"particle_swarm",
                        "swarm_size":300,
                        "iteration_limit":6000}}
                
    def extractModelParam(self):
        copasi_filename = self.genPathCopasi("extractor")
        self.recentModel = model.loada(self.antString, copasi_filename)
        return self.recentModel.parameters.copy().squeeze().to_dict()
            
    def genPathCopasi(self,nameBase,suffix=".cps"):
        """generates a copasi path that isn't being used yet
        
        we want to avoid acidently useing the same copasi file for difrent
        calculations as it seems they don't always get overwriten. This
        function generates a path to an unused file name in the classes
        run directory by trying difrent sufixes.
        
        Args:
           nameBase (str):  base to use in trying to creat file name.
           
        Kwargs:
           suffix (str): optional string for file extention incase we ever
               want to generate a path to non copasi file someday. defaults
               to ".cps"
               
        Returns:
           str: path to unused filename in objects run directory.
        """
        i=0
        nameFree=False
        while not nameFree:
            copasi_filename = os.path.join(self.run_dir,nameBase+
                                           str(i)+suffix)
            nameFree = not os.path.exists(copasi_filename)
            i=i+1
        return copasi_filename
    
    def clearRunDirectory(self):
        """empties the objects run directory
        
        empties the objects run directory of all '.cps', '.txt', '.sbml',
        '.csv' extention files and any empty directories. It seems
        pycotools creats a lot of files in paramiter estimations it has a
        habit of reusing. Clearing the run directory is one way to prevent
        this.
        """
        for root, dirs, files in os.walk(self.run_dir, topdown=False):
            for name in files:
                if name.lower().endswith(('.cps', '.txt', '.sbml', '.csv')):
                    os.remove(os.path.join(root, name))
            for name in dirs:
                if len(os.listdir(os.path.join(root, name)))==0:
                    os.rmdir(os.path.join(root, name))
            
    def genPrefixAntString(self,estimatedVar,prefix="_"):
        """Creates a new antimony string with added prefixs
        
        creats a new antimony string (stored seperatly from the refrence
        string) where certain variables in the string have had a prefix
        added. Pycotools uses prefixes to determine which variables should
        be estimated in a parameter estimation.
        
        Args:
           estimatedVar (list of str):  the variables in the antimony
               string (and hence the model) that should be prefixed.
               
        Kwargs:
           prefix (str): optionaly the prefix to add
        """
        self.prefixAntString = self.antString
        for name in estimatedVar:
            self.prefixAntString = replaceVariable(self.prefixAntString,
                                                   name,prefix+name)
    
    def genRefCopasiFile(self, filePath = None):
        if filePath is None:
            copasi_filename = self.genPathCopasi("refFile")
        else:
            copasi_filename = filePath
        self.recentModel = model.loada(self.antString, copasi_filename)
        
    def runSteadyStateFinder(self,params=None):
        r = te.loada(self.antString)
        outValues = r.getFloatingSpeciesIds()
        boundValues = r.getBoundarySpeciesIds()
        outValues.extend(boundValues)
        r.conservedMoietyAnalysis = True
        r.setSteadyStateSolver('nleq1')
        if params is not None:
            for key, value in params.items():
                temp = r[key]
                r[key] = value
                print(key,": ",temp," -> ",r[key])
        
        r.steadyState()
        rState = {key:r[key] for key in outValues}
        return rState
    
    # if you reuse the same PEName twise the results get over writen by
    # the first
    def runParamiterEstimation(self,expDataFP,PEName=None,
                               copyNum=1,estimateIC=True,estimatedVar=None,
                               prefix='_',rocket=False,upperParamBound=None,
                               lowerParamBound=None,
                               method="particle_swarm_default",
                               overrideParam=None,
                               indepToAdd=None, endTime=None):
        if PEName is None:
            copasi_filename = self.genPathCopasi("paramiterEstimation")
        else:
            copasi_filename = os.path.join(self.run_dir, PEName)
        
        if estimatedVar is not None:
            self.genPrefixAntString(estimatedVar)
            self.recentModel = model.loada(self.prefixAntString,
                                           copasi_filename)
            colRenameDict = {i:(prefix+i) for i in estimatedVar}
            expDataFP = renameCSVColumns(expDataFP,self.run_dir,
                                         colRenameDict,
                                         indepToAdd=indepToAdd)
        else:
            self.recentModel = model.loada(self.antString, copasi_filename)
        if overrideParam is not None:
            model.InsertParameters(self.recentModel,
                                   parameter_dict=overrideParam,
                                   inplace=True)
        
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
            if upperParamBound is not None:
                context.set('upper_bound', upperParamBound)
            if lowerParamBound is not None:
                context.set('lower_bound', lowerParamBound)
            for key, value in self.methodDict[method].items():
                context.set(key, value)
            context.set('run_mode', runMode) 
            context.set('pe_number', 1) 
            context.set('copy_number', copyNum) 
            config = context.get_config()
        
        self.recentPE = tasks.ParameterEstimation(config)
        doneNum=0
        lastLoop=0
        logging.disable(logging.WARNING)
        while copyNum>doneNum:
            try:
                parse_object = viz.Parse(self.recentPE)
                if isinstance(parse_object[list(parse_object)[0]],
                                          pd.DataFrame):
                    doneNum = parse_object[list(parse_object)[0]].shape[0]
            except:
                time.sleep(60)
            if rocket and lastLoop<doneNum:
                print(str(doneNum) + ' of ' + str(copyNum) + ' done')
                lastLoop = doneNum
            if endTime is not None and rocket:
                if endTime<=time.time():
                    break
        
        logging.disable(logging.NOTSET)
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
    
    def preProcessParamEnsam(self,Params):
        copasi_filename = self.genPathCopasi("adjuster")
        self.recentModel = model.loada(self.antString, copasi_filename)
        adjust = {}
        for colName in list(Params.columns):
            baseVal = self.recentModel.get('metabolite', colName, by='name')
            baseVal = baseVal.to_df()
            adjust[colName] = baseVal["Value"]["concentration"]
        return Params.fillna(value=adjust)
    
    def TCendState(self, timeCourse, variables = None):
        if isinstance(timeCourse,list):
            myTC=timeCourse
        else:
            myTC=[timeCourse]
        myTC = [i.tail(1).squeeze().to_dict() for i in myTC]
        myTC = pd.DataFrame(myTC)
        myTC = myTC.drop(columns=["Time"])
        if variables == "metabolites":
            copasi_filename = self.genPathCopasi("getter")
            self.recentModel = model.loada(self.antString, copasi_filename)
            myVar = [i.name for i in self.recentModel.metabolites]
        if variables is not None:
            myTC = myTC[myVar]
        return myTC
    
    def runTimeCourse(self,duration,stepSize=0.01,intervals=100,
                      TCName=None,adjustParams=None,subSet=None,
                      rocket=False):
        if adjustParams is not None:
            if subSet is None:
                subSet = range(len(adjustParams.index))
            results = []
            if not isinstance(duration,list):
                duration = [duration for i in subSet]
            if len(duration)!=len(subSet):
                return False
            if rocket:
                timeCourses = []
                for setIndex, myDur in zip(subSet,duration):
                    copasi_filename = self.genPathCopasi("timeCourse")
                    self.recentModel = model.loada(self.antString,
                                                   copasi_filename)
                    model.InsertParameters(self.recentModel,
                                           df=adjustParams,
                                           index=setIndex,inplace=True)
                    self.recentTimeCourse = tasks.TimeCourse(
                            self.recentModel,end=myDur,
                            step_size=stepSize,intervals=intervals,
                            run=False)
                    timeCourses.append(self.recentTimeCourse)
                    myScriptName = self.genPathCopasi("TCSlurmScript",
                                                      suffix = ".sh")
                    shellString = ("#!/bin/bash\nCopasiSE "+
                                   copasi_filename)
                    f = open(myScriptName, 'w')
                    f.write(shellString)
                    f.close()
                    os.system("sbatch "+myScriptName)
                for theTimeCourse in timeCourses:
                    sucsessful=False
                    logging.disable(logging.WARNING)
                    while not sucsessful:
                        try:
                            parse_object = viz.Parse(theTimeCourse)
                            results.append(parse_object.data.copy())
                            sucsessful = True
                        except:
                            time.sleep(60)
                    logging.disable(logging.NOTSET)
                return results
            else:
                if TCName is None:
                    copasi_filename = self.genPathCopasi("timeCourse")
                else:
                    copasi_filename = os.path.join(self.run_dir, TCName)
                self.recentModel = model.loada(self.antString,
                                               copasi_filename)
                for setIndex, myDur in zip(subSet,duration):
                    model.InsertParameters(self.recentModel,
                                           df=adjustParams,
                                           index=setIndex,inplace=True)
                    self.recentTimeCourse = tasks.TimeCourse(
                            self.recentModel,end=myDur,
                            step_size=stepSize,intervals=intervals)
                    results.append(viz.Parse(
                            self.recentTimeCourse).data.copy())
            return results
        else:
            if TCName is None:
                copasi_filename = self.genPathCopasi("timeCourse")
            else:
                copasi_filename = os.path.join(self.run_dir, TCName)
                
            self.recentModel = model.loada(self.antString, copasi_filename)
            self.recentTimeCourse = tasks.TimeCourse(self.recentModel,
                                                     end=duration,
                                                     step_size=stepSize,
                                                     intervals=intervals)
            return viz.Parse(self.recentTimeCourse).data.copy()