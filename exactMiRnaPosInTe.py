#!/usr/bin/env python3

import sys
from Bio import SeqIO

# find exact possions of all miRNA within TEs and output summary tables 
# specific for given TE family

miRnasProtDomBed = sys.argv[1]
allProtDomTesGff = sys.argv[2]

def genTeDict(superFamList):
    teDict = {}
    for frag in superFamList:
        teDict[frag] = [0, 0]
    return teDict

def getTeHeaderList(miRnasProtDomBed):
    teHeadList = []
    with open(miRnasProtDomBed) as bed:
        for rec in bed:
            teHeadList.append(rec.split("\t")[0])
    teHeadListUniq = list(set(teHeadList))
    return teHeadListUniq

def fillTeDict(teBasicDict, te, allProtDomTesGff):
    # !!! Tady sem skoncil !!!

def getMiRnaPos(miRnasProtDomBed, allProtDomTesGff):
    # upload coordinates of resp. TE fragments into dictionary
    copiaList = ['ltr left', 'pbs', 'GAG', 'PROT', 'INT', 'RT', 'RH', 'ppt', 'ltr right']
    gypsyChrList = ['ltr left', 'pbs', 'GAG', 'PROT', 'RT', 'RH', 'INT', 'CHD', 'ppt', 'ltr right']
    gypsyNonChrList = ['ltr left', 'pbs', 'GAG', 'PROT', 'RT', 'RH', 'INT', 'ppt', 'ltr right']
    # copiaListX = ['ltr left', 'ltr left-pbs', 'pbs', 'pbs-GAG', 'GAG', 'GAG-PROT', 'PROT', 'PROT-INT', 'INT', 
    # 'INT-RT', 'RT', 'RT-RH', 'RH', 'RH-ppt', 'ppt', 'ppt-ltr right', 'ltr right']
    # gypsyChrListX = ['ltr left', 'ltr left-pbs', 'pbs', 'pbs-GAG', 'GAG', 'GAG-PROT', 'PROT', 'PROT-RT', 'RT', 
    # 'RT-RH', 'RH', 'RH-INT', 'INT', 'INT-CHD', 'CHD', 'CHD-ppt', 'ppt', 'ppt-ltr right', 'ltr right']
    # gypsyNonChrListX = ['ltr left', 'ltr left-pbs', 'pbs', 'pbs-GAG', 'GAG', 'GAG-PROT', 'PROT', 'PROT-RT', 'RT', 
    # 'RT-RH', 'RH', 'RH-INT', 'INT', 'INT-ppt', 'ppt', 'ppt-ltr right', 'ltr right']
    copiaDict = genTeDict(copiaList)
    gypsyChrDict = genTeDict(gypsyChrList)
    gypsyNonChrDict = genTeDict(gypsyNonChrList)

    # here I have to define which families belong to which superFams:
    copiaFams = ["copia_Ale","copia_Alesia","copia_Ivana","copia_Lyco","copia_Osser","copia_SIRE","copia_TAR",
    "copia_Tork","copia_Ty1-outgroup"]
    gypsyChrFams = ["gypsy_Chlamyvir","gypsy_chromo-outgroup","gypsy_chromo-unclass","gypsy_CRM",
    "gypsy_Galadriel","gypsy_Reina","gypsy_Tcn1","gypsy_Tekay"]
    gypsyNonChrFams = ["gypsy_non-chromo-outgroup","gypsy_Ogre","gypsy_Retand","gypsy_Selgy","gypsy_TatI","gypsy_TatII",
    "gypsy_TatIII","gypsy_Athila","gypsy_Phygy","gypsy_Selgi"]
    famList = [copiaFams, gypsyChrFams, gypsyNonChrFams]
    teHeadList = getTeHeaderList(miRnasProtDomBed)
    
    for te in teHeadList:
        # zjistiti fam a vytvorit dictionary
        fam = te.split("|")[3]
        teBasicDict = {}    # bez prostoru mezi fragmenty
        for famGroup in famList:
            if fam in famGroup:
                teBasicDict = genTeDict(famGroup)
                teBasicDict = fillTeDict(teBasicDict, te, allProtDomTesGff)
        teFullDict = {}     # se zaclenenim prostoru mezi fragmenty a jejich koordinat

    with open(miRnasProtDomBed) as bed:
        for rec in bed:
            header = rec.split("\t")[0]
            with open(allProtDomTesGff) as gff:
                for l in gff:
                    if header in l:
                        print(header,l)

getMiRnaPos(miRnasProtDomBed, allProtDomTesGff)
