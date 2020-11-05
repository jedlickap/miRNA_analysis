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
        print("XXX")
        antiFragList = list(teDict.keys())
        coordList = []
        for k in teDict.keys():
            coordList.append(teDict[k][1])
            coordList.append(teDict[k][0])
        print("coordList", coordList)
        teNewDict = {}
        for i in range(len(antiFragList)):
            lastValue = 0
            if i == 0:
                teNewDict[antiFragList[i]] = []
                teNewDict[antiFragList[i]].append(1)
                teNewDict[antiFragList[i]].append((coordList[0] - coordList[1]) + 1)
                lastValue = teNewDict[antiFragList[i]][-1]
                del coordList[0]
            elif len(coordList) != 3:
                teNewDict[antiFragList[i]] = []
                teNewDict[antiFragList[i]].append((coordList[0] - coordList[1]) + lastValue)
                teNewDict[antiFragList[i]].append((coordList[1] - coordList[2]) + teNewDict[antiFragList[i]][-1])
                lastValue = teNewDict[antiFragList[i]][-1]
                del coordList[0]
            else:
                teNewDict[antiFragList[i]] = []
                teNewDict[antiFragList[i]].append((coordList[0] - coordList[1]) + lastValue)
                teNewDict[antiFragList[i]].append(1)  
        teDict = teNewDict      
    return teDict

def getMiRnaPos(miRnasProtDomBed, allProtDomTesGff):
    # copiaListX = ['ltr left', 'ltr left-pbs', 'pbs', 'pbs-GAG', 'GAG', 'GAG-PROT', 'PROT', 'PROT-INT', 'INT', 
    # 'INT-RT', 'RT', 'RT-RH', 'RH', 'RH-ppt', 'ppt', 'ppt-ltr right', 'ltr right']
    # gypsyChrListX = ['ltr left', 'ltr left-pbs', 'pbs', 'pbs-GAG', 'GAG', 'GAG-PROT', 'PROT', 'PROT-RT', 'RT', 
    # 'RT-RH', 'RH', 'RH-INT', 'INT', 'INT-CHD', 'CHD', 'CHD-ppt', 'ppt', 'ppt-ltr right', 'ltr right']
    # gypsyNonChrListX = ['ltr left', 'ltr left-pbs', 'pbs', 'pbs-GAG', 'GAG', 'GAG-PROT', 'PROT', 'PROT-RT', 'RT', 
    # 'RT-RH', 'RH', 'RH-INT', 'INT', 'INT-ppt', 'ppt', 'ppt-ltr right', 'ltr right']
    
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
    
    for te in teHeadList:
        # zjistiti fam a vytvorit dictionary
        fam = te.split("|")[3]        
        for superFam in sorted(famFragsDict.keys()):
            if fam in famFragsDict[superFam][0]:
                teBasicDict = fillTeDict(famFragsDict[superFam][1], te, allProtDomTesGff)
                for k in teBasicDict.keys():
                    print(k, teBasicDict[k])

    #     teFullDict = {}     # se zaclenenim prostoru mezi fragmenty a jejich koordinat

    # with open(miRnasProtDomBed) as bed:
    #     for rec in bed:
    #         header = rec.split("\t")[0]
    #         with open(allProtDomTesGff) as gff:
    #             for l in gff:
    #                 if header in l:
    #                     print(header,l)

getMiRnaPos(miRnasProtDomBed, allProtDomTesGff)
