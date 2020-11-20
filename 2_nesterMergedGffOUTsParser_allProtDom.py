#!/usr/bin/env python3

import sys
import re
from operator import itemgetter
from Bio import SeqIO
import subprocess
import os

miRnaBed = sys.argv[1]
gff = sys.argv[2]
tesFa = sys.argv[3] # TEs fasta with '>Gmax|Chr01|TE_827|gypsy_CRM|noTSD' headers.
miRnaFa = sys.argv[4]
preMiRnaFa = sys.argv[5]

"""
long_terminal_repeat [ltr]: left, right
nested_repeat [nr]
polypeptide_conserved_region [pcr]: aRH,CHD,CHDCR,GAG,INT,PROT,RH,RT.
primer_binding_site [pbs]
repeat_fragment [rf]
RR_tract [ppt]
target_site_duplication [tsd]: left, right
### all features will be in list["start", "end", "8th column from GFF"]
"""
def mkdirForOut():
    outFolder = 'tesWithAllProtDom_outputs'
    if not os.path.exists(outFolder):
        os.makedirs(outFolder, mode=0o777)        
    return outFolder + '/'

def fillFeatureList(l):
    l = l.rstrip()
    fList = []
    fList.append(l.split("\t")[3])
    fList.append(l.split("\t")[4])
    fList.append(l.split("\t")[8])
    return fList

def fillInDict(gff):
    inDict = {}
    reg = r'=(TE_BASE \d+);'
    chrom = ""
    te = ""
    with open(gff) as f:        
        for l in f:
            if "##" not in l:                    
                if "nested_repeat" in l:
                    chrom = str(l.split("\t")[0])
                    te = re.search(reg, l).group(1).replace("TE_BASE ", "")
                    if chrom not in inDict.keys():
                        inDict[chrom] = {}
                        inDict[chrom][te] = {}
                        inDict[chrom][te]["nested_repeat"] = []
                        inDict[chrom][te]["nested_repeat"] = fillFeatureList(l)
                    else:
                        inDict[chrom][te] = {}
                        inDict[chrom][te]["nested_repeat"] = []
                        inDict[chrom][te]["nested_repeat"] = fillFeatureList(l)                    
                if chrom in l and te + ";" in l and "nested_repeat" not in l:
                    if "repeat_fragment" in l:
                        if "repeat_fragment" not in inDict[chrom][te].keys():
                            inDict[chrom][te]["repeat_fragment"] = []
                            inDict[chrom][te]["repeat_fragment"].append(fillFeatureList(l))
                        else:
                            inDict[chrom][te]["repeat_fragment"].append(fillFeatureList(l))
                    if "polypeptide_conserved_region" in l:
                        if "polypeptide_conserved_region" not in inDict[chrom][te].keys():
                            inDict[chrom][te]["polypeptide_conserved_region"] = []
                            inDict[chrom][te]["polypeptide_conserved_region"].append(fillFeatureList(l))
                        else:
                            inDict[chrom][te]["polypeptide_conserved_region"].append(fillFeatureList(l))
                    if "primer_binding_site" in l:
                        if "primer_binding_site" not in inDict[chrom][te].keys():
                            inDict[chrom][te]["primer_binding_site"] = ""
                            inDict[chrom][te]["primer_binding_site"] = fillFeatureList(l)
                        else:
                            inDict[chrom][te]["primer_binding_site"] = fillFeatureList(l)
                    if "RR_tract" in l:
                        if "RR_tract" not in inDict[chrom][te].keys():
                            inDict[chrom][te]["RR_tract"] = ""
                            inDict[chrom][te]["RR_tract"] = fillFeatureList(l)
                        else:
                            inDict[chrom][te]["RR_tract"] = fillFeatureList(l)
                    if "long_terminal_repeat" in l:
                        if "long_terminal_repeat" not in inDict[chrom][te].keys():
                            inDict[chrom][te]["long_terminal_repeat"] = []
                            inDict[chrom][te]["long_terminal_repeat"].append(fillFeatureList(l))
                        else:
                            inDict[chrom][te]["long_terminal_repeat"].append(fillFeatureList(l))
                    if "target_site_duplication" in l:
                        if "target_site_duplication" not in inDict[chrom][te].keys():
                            inDict[chrom][te]["target_site_duplication"] = []
                            inDict[chrom][te]["target_site_duplication"].append(fillFeatureList(l))
                        else:
                            inDict[chrom][te]["target_site_duplication"].append(fillFeatureList(l))
    return inDict

featureList = ["nested_repeat", "repeat_fragment", "polypeptide_conserved_region", "primer_binding_site", "RR_tract", "long_terminal_repeat", "target_site_duplication"]

def getAllDomsDict(gff, inDict, miRnaBed, outFolder):
    cnt = 0
    allProtDomTesBed_name = miRnaBed.replace(".bed","_allProtDomTes.bed")
    bedOutList = []
    allDomsDict = {}
    with open(miRnaBed) as headKeys:
        next(headKeys)
        for rec in headKeys:
            recTrue = False
            ch, teID = rec.split("\t")[0].split("|")[1], rec.split("\t")[0].split("|")[2].split("_")[1]
            if 'polypeptide_conserved_region' in inDict[ch][teID].keys():
                if inDict[ch][teID]['polypeptide_conserved_region']:
                    domDict = {'RH;color':True,'name=CHD':True,'name=GAG':True,'name=INT':True,'name=PROT':True,'name=RH':True}
                    for dom in sorted(domDict.keys()):
                        for domRec in inDict[ch][teID]['polypeptide_conserved_region']:
                            if dom in "\t".join(domRec):
                                domDict[dom] = True
                                break
                            else:
                                domDict[dom] = False
                    if False not in list(domDict.values()):
                        #print(ch, teID, inDict[ch][teID])
                        recTrue = True
                    domDict['name=CHD'] = True
                    if False not in list(domDict.values()):
                        #print(ch, teID, inDict[ch][teID])
                        recTrue = True
                    if recTrue:
                        cnt +=1
                if recTrue:
                    bedOutList.append(rec)
                    if ch not in allDomsDict.keys():
                        allDomsDict[ch] = {}
                        if teID not in allDomsDict[ch].keys():
                            allDomsDict[ch][teID] = {}
                            allDomsDict[ch][teID] = inDict[ch][teID]
                    else:
                        if teID not in allDomsDict[ch].keys():
                            allDomsDict[ch][teID] = {}
                            allDomsDict[ch][teID] = inDict[ch][teID]
    currDir = os.getcwd()
    os.chdir(outFolder)
    with open(allProtDomTesBed_name,"w") as bedOut:
        for l in bedOutList:
            bedOut.write(l)
    os.chmod(allProtDomTesBed_name, 0o777)                        
    teCnt = 0
    os.chdir(currDir)
    for ch in allDomsDict.keys():
        teCnt += len(allDomsDict[ch].keys())        
    print("Number of TEs with all protDoms in {}: {} ({} miRNAs)".format(gff, teCnt, cnt))
    return allDomsDict, allProtDomTesBed_name

def getBase(gff, ch, teIDwhole):
    baseNr = 0
    teLen = 0
    gffPat = ch + "\tfeature\tnested_repeat\t"
    with open(gff) as gffIn:
        for l in gffIn:
            if gffPat in l and '=' + teIDwhole + ';' in l:
                baseNr = int(l.split('\t')[3]) - 1
                teLen = int(l.split('\t')[4]) - baseNr
    return baseNr, teLen

def printOutFasta(sourseFa, faName, headListUniq):
    with open(faName, "w") as outFa:    
        for head in headListUniq:
            for fa in SeqIO.parse(sourseFa, "fasta"):
                if head == fa.id or head in fa.id:
                    outFa.write(">" + fa.id + '\n')
                    outFa.write(str(fa.seq) + "\n")
    os.chmod(faName, 0o777)

def getAllDomsGffFa(allProtDomTesBed_name, gff, tesFa, miRnaFa, preMiRnaFa, outFolder):
    gffOutName = outFolder + allProtDomTesBed_name.replace(".bed",".gff")
    faOutName = outFolder + allProtDomTesBed_name.replace(".bed",".fa")
    miRnaFaOutName = outFolder + miRnaFa.replace(".fa","_allProtDomTes.fa")
    preMiRnaFaOutName =outFolder + preMiRnaFa.replace(".fa","_allProtDomTes.fa")
    allProtDomTesBed_name = outFolder + allProtDomTesBed_name
    headList = []
    with open(allProtDomTesBed_name) as allProtDomTesBed:
        for rec in allProtDomTesBed:
            headList.append(rec.split("\t")[0])
    headListUniq = list(set(headList))
    # generate allProtDomTesFa  
    printOutFasta(tesFa, faOutName, headListUniq)
    printOutFasta(miRnaFa, miRnaFaOutName, headListUniq)
    printOutFasta(preMiRnaFa, preMiRnaFaOutName, headListUniq)
    # generate allProtDomTesBed
    with open(gffOutName, "w") as out:
        for head in headListUniq:
            ch, teID = head.split("\t")[0].split("|")[1], head.split("\t")[0].split("|")[2].split("_")[1]
            teIDwhole = "TE_BASE " + str(teID)
            baseNr, teLen = getBase(gff, ch, teIDwhole)
            out.write("##gff-version 3\n")
            out.write("##sequence-region {} 1 {}\n".format(head, teLen))
            with open(gff) as gffIn:
                for l in gffIn:
                    if ch in l and '=' + teIDwhole + ';' in l and 'target_site_duplication' not in l:
                        lineList = l.split('\t')
                        lineList[3] = str(int(lineList[3]) - baseNr)
                        lineList[4] = str(int(lineList[4]) - baseNr)
                        out.write(head + "\t" + "\t".join(lineList[1:]))
    os.chmod(gffOutName, 0o777)
    return gffOutName

# create output folder
outFolder = mkdirForOut()      
# get Gff into dictionary
inDict = fillInDict(gff)
# print out bed file for miRNA in TEs with all prot. domain only filtered gff dictionary 
allDomsDict, allProtDomTesBed_name = getAllDomsDict(gff, inDict, miRnaBed, outFolder)
# print gff and fasta file filtered for TEs with miRNAs, miRNAs, premiRNAs all with all prot. domain only filtered gff dictionary 
gffOutName = getAllDomsGffFa(allProtDomTesBed_name, gff,tesFa, miRnaFa, preMiRnaFa, outFolder)
