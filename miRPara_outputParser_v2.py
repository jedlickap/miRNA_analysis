#!/usr/bin/env python3

######
# this program should parse miRPara '*_level_1.out' file in order to find:
# 1. best fitted predictions of miRNAs
# 2. recalculate their respective position within genomic sequences

import sys
import re

miRParaOut = sys.argv[1]
headerKeys = sys.argv[2]

# upload miRPara output into dictionary:
def getMirParaDict(miRParaOut):
    miRParaDict = {}
    with open(miRParaOut) as inPut:
        for l in inPut:
            if not re.search(r'^#', l):
                recList = l.rstrip().split("\t")
                if recList[0] not in miRParaDict.keys():
                    miRParaDict[recList[0]] = {}
                    miRParaDict[recList[0]]['5P'] = []
                    miRParaDict[recList[0]]['3P'] = []
                    if "\t5P\t" in l:
                        miRParaDict[recList[0]]['5P'].append("\t".join(recList[1:]))
                    if "\t3P\t" in l:
                        miRParaDict[recList[0]]['3P'].append("\t".join(recList[1:]))
                else:
                    if "\t5P\t" in l:
                        miRParaDict[recList[0]]['5P'].append("\t".join(recList[1:]))
                    if "\t3P\t" in l:
                        miRParaDict[recList[0]]['3P'].append("\t".join(recList[1:]))
    return miRParaDict

def findHighestSvmProb(recList):
    svmProbList = []
    for r in recList:
        svmProbList.append(float(r.split("\t")[4]))
    maxIndex = svmProbList.index(max(svmProbList))
    return str(recList[maxIndex])

def filterCandidates(inDict):
    bestMiRnasDict = {}
    for te in sorted(inDict.keys()):
        bestMiRnasDict[te] = {}
        for primeEnd in sorted(inDict[te].keys()):
            if inDict[te][primeEnd]:
                bestMiRnasDict[te][primeEnd] = findHighestSvmProb(inDict[te][primeEnd])
            else:
                bestMiRnasDict[te][primeEnd] = None
    return bestMiRnasDict

def getHeader(teKey, headerKeys):
    te = teKey.split(":")[0]
    header = None
    with open(headerKeys) as hk:
        for l in hk:
            if "\t" + te + "\n" in l:
                header = l.split("\t")[0]
    return header

def getStEnDescr(oneMiRnaItem, headerKeys):
    recList = oneMiRnaItem.split("\t")
    teSeList = recList[1].split(":")
    seList = teSeList[1].split("-")
    seList[0] = int(seList[0]) + (int(teSeList[2].split("_")[0]) - 2)
    seList[1] = int(seList[0]) + int(teSeList[2].split("_")[1])
    faHeader = getHeader(teSeList[0], headerKeys)
    return "{}\t{}\t{}\t{}\n".format(faHeader, str(seList[0]), str(seList[1]), ";".join(recList[1:]))

def addHeaderOutput(bestMiRnasDict, headerKeys, miRParaOut): 
    # print out file of filtered miRNAs with the highest SVMprob
    # TODO print out gff/bed file for visualization and intersect determination!!!
    outName = miRParaOut.replace(".out","_filtered.out")
    bedName = miRParaOut.replace(".out","_filtered.bed")
    with open(outName, "w") as out:
        with open(bedName, "w") as out1:
            out1.write("track name=\"miRNA in {}\" color=255,0,0\n".format(miRParaOut))
            for k in sorted(bestMiRnasDict.keys()):
                for p in sorted(bestMiRnasDict[k].keys()):
                    if bestMiRnasDict[k][p]:
                        out.write("{}\t{}\t{}\n".format(getHeader(k, headerKeys), k, bestMiRnasDict[k][p]))
                        out1.write(getStEnDescr(bestMiRnasDict[k][p], headerKeys))

def main(miRParaOut, headerKeys):
    inDict = getMirParaDict(miRParaOut) # upload miRPara output
    bestMiRnasDict = filterCandidates(inDict) # filter miRNA with the best SVM_score
    cnt = 0
    for k in bestMiRnasDict.keys():
        for p in bestMiRnasDict[k]:
            if bestMiRnasDict[k][p]:
                cnt += 1
    print("Nr. of filtered miRNAs: {}".format(cnt))
    addHeaderOutput(bestMiRnasDict, headerKeys, miRParaOut) # print results to the text file with original headers in the 1st column: 'Acom|LG01|TE_132|gypsy_Tekay|noTSD'

if __name__ == '__main__':
    main(miRParaOut, headerKeys)
