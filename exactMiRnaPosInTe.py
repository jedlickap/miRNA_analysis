#!/usr/bin/env python3

import sys
from Bio import SeqIO

# find exact possions of all miRNA within TEs and output summary tables 
# specific for given TE family

miRnasProtDomBed = sys.argv[1]
allProtDomTesGff = sys.argv[2]

def getTeHeaderList(miRnasProtDomBed):
    teHeadList = []
    with open(miRnasProtDomBed) as bed:
        for rec in bed:
            teHeadList.append(rec.split("\t")[0])
    teHeadListUniq = list(set(teHeadList))
    return teHeadListUniq

def antisenseRecalc(teDict):
    antiFragList = list(teDict.keys())
    antiFragList.reverse()
    coordList = []
    for k in teDict.keys():
        coordList.append(teDict[k][1])
        coordList.append(teDict[k][0])
    teNewDict = {}
    
    for i in range(len(antiFragList)):
        if i == 0:
            teNewDict[antiFragList[i]] = []
            teNewDict[antiFragList[i]].append(1)
            teNewDict[antiFragList[i]].append(1 + (coordList[0] - coordList[1]))
            lastValue = teNewDict[antiFragList[i]][1]
            del coordList[0]
        elif len(coordList) != 3:
            teNewDict[antiFragList[i]] = []
            teNewDict[antiFragList[i]].append((coordList[0] - coordList[1]) + lastValue)
            teNewDict[antiFragList[i]].append((coordList[1] - coordList[2]) + teNewDict[antiFragList[i]][0])
            lastValue = teNewDict[antiFragList[i]][1]
            del coordList[0]
            del coordList[0]
        else:
            teNewDict[antiFragList[i]] = []
            teNewDict[antiFragList[i]].append((coordList[0] - coordList[1]) + lastValue)
            teNewDict[antiFragList[i]].append((coordList[1] - coordList[2]) + teNewDict[antiFragList[i]][0])  
    teDict = teNewDict
    return teDict

def fillTeDict(fragList, te, allProtDomTesGff):
    teDict = {}
    for i in range(len(fragList)):
        pat = ";name=" + fragList[i] + ";"
        with open(allProtDomTesGff) as gff:
            for l in gff:
                if te in l and pat in l:
                    teDict[fragList[i]] = []
                    teDict[fragList[i]].append(int(l.split("\t")[3]))
                    teDict[fragList[i]].append(int(l.split("\t")[4]))   
    # recalculate coordinates if antisense TE = True
    if teDict['ltr left'][0] > teDict['ltr right'][0]:
        teDict = antisenseRecalc(teDict)              
    return teDict

def posFinder(start, end, teBasicDict):
    found = False
    outList = [] # fragName, fragS, fragE
    fragList = list(teBasicDict.keys())
    # check inside regular fragments
    for k in teBasicDict.keys():
        if not found:
            if start >= teBasicDict[k][0] and end <= teBasicDict[k][1]:
                found = True  
                outList = [k, teBasicDict[k][0], teBasicDict[k][1]]
                break
    # check within the space between two adjacent fragments
    # and
    # check which two parts does TE miRNA intersect
    if not found:
        for i in range(len(fragList)):
            if fragList[i] != fragList[-1]:
                if teBasicDict[fragList[i]][1] <= start and end <= teBasicDict[fragList[i+1]][0]:                
                    outList = [fragList[i] + "-" + fragList[i+1], teBasicDict[fragList[i]][1], 
                    teBasicDict[fragList[i+1]][0]]
                    found = True 
                    break
                elif teBasicDict[fragList[i]][1] > start and (teBasicDict[fragList[i]][1] < end < teBasicDict[fragList[i+1]][0]):
                    outList = [fragList[i] + ":" + fragList[i] + "-" + fragList[i+1], 
                    teBasicDict[fragList[i]][0], teBasicDict[fragList[i+1]][0]]
                    found = True
                    break
                elif (teBasicDict[fragList[i]][1] < start < teBasicDict[fragList[i+1]][0]) and end > teBasicDict[fragList[i+1]][0]:
                    outList = [fragList[i] + "-" + fragList[i+1] + ":" + fragList[i+1], 
                    teBasicDict[fragList[i]][1], teBasicDict[fragList[i+1]][1]]
                    found = True
                    break                
            else:
                if start < teBasicDict[fragList[i]][0] and end < teBasicDict[fragList[i]][1]:
                    outList = [fragList[i-1] + "-" + fragList[i] + ":" + fragList[i], 
                    teBasicDict[fragList[i-1]][1], teBasicDict[fragList[i]][1]]
                    found = True
                    break
    return "\t".join([str(elem) for elem in outList])   # join([str(elem) for elem in s]) 

def miRnaCoordInTe(te, teBasicDict, miRnasProtDomBed):
    teMiRnasOutList = []
    miRnasList = []
    with open(miRnasProtDomBed) as bed:
        for rec in bed:
            if te in rec:
                miRnasList.append(rec)
    for miRna in miRnasList:
        start, end = int(miRna.split("\t")[1]), int(miRna.split("\t")[2])
        teMiRnasOutList.append(miRna.split("\t")[0].replace("|", "\t") + "\t" + posFinder(start, end, teBasicDict) + 
        "\t" + str(start) + "\t" + str(end) + "\t" + miRna.split("\t")[3].rstrip())                                        # fragName, fragLen, miRnaRelPos
    return teMiRnasOutList        

def getMiRnaPos(miRnasProtDomBed, allProtDomTesGff):
    # here I have to define which families belong to which superFams:
    famFragsDict = {
        'copia':[
            ["copia_Ale","copia_Alesia","copia_Ivana","copia_Lyco","copia_Osser","copia_SIRE","copia_TAR","copia_Tork","copia_Ty1-outgroup"], 
            ['ltr left', 'pbs', 'GAG', 'PROT', 'INT', 'RT', 'RH', 'ppt', 'ltr right']], 
        'gypsyChr': [
            ["gypsy_Chlamyvir","gypsy_chromo-outgroup","gypsy_chromo-unclass","gypsy_CRM","gypsy_Galadriel","gypsy_Reina","gypsy_Tcn1","gypsy_Tekay"], 
            ['ltr left', 'pbs', 'GAG', 'PROT', 'RT', 'RH', 'INT', 'CHD', 'ppt', 'ltr right']],
        'gypsyNonChr': [
            ["gypsy_non-chromo-outgroup","gypsy_Ogre","gypsy_Retand","gypsy_Selgy","gypsy_TatI","gypsy_TatII","gypsy_TatIII","gypsy_Athila","gypsy_Phygy",
            "gypsy_Selgi"], 
            ['ltr left', 'pbs', 'GAG', 'PROT', 'RT', 'RH', 'INT', 'ppt', 'ltr right']]
            }
    teHeadList = getTeHeaderList(miRnasProtDomBed)
    outputList = []

    # TODO: include output of all TEs fragments lengths and lengths of distances between them!

    for te in teHeadList:
        # which TE family and add specific fragments order
        fam = te.split("|")[3]        
        for superFam in sorted(famFragsDict.keys()):
            if fam in famFragsDict[superFam][0]:
                # generate dictionary with all coordinates of TE fragments
                teBasicDict = fillTeDict(famFragsDict[superFam][1], te, allProtDomTesGff)
                # here to add line prep for TE frag lengths!
                outputList.append(miRnaCoordInTe(te, teBasicDict, miRnasProtDomBed))
    
    # get tsv tab le output:
    spec = miRnasProtDomBed.split("_")[0]
    outName = spec + "_miRNA_accuratePosInTe.tsv"
    with open(outName, "w") as out:
        for l1 in outputList:
            for l2 in l1:
                out.write(l2 + "\n")
                        
getMiRnaPos(miRnasProtDomBed, allProtDomTesGff)
