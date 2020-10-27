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
    # TODO print out fasta files with both the pre-miRNAs and miRNAs
    outName = miRParaOut.replace(".out","_filtered.out")
    bedName = miRParaOut.replace(".out","_filtered.bed")
    preMiRnaFaName = miRParaOut.replace(".out","_preMiRna.fa")
    miRnaFaName = miRParaOut.replace(".out","_miRna.fa")

    with open(outName, "w") as out:
        with open(bedName, "w") as out1:
            with open(preMiRnaFaName, "w") as out2:
                with open(miRnaFaName, "w") as out3:
                    out1.write("track name=\"miRNA in {}\" color=255,0,0\n".format(miRParaOut))
                    for k in sorted(bestMiRnasDict.keys()):
                        for p in sorted(bestMiRnasDict[k].keys()):
                            if bestMiRnasDict[k][p]:
                                out.write("{}\t{}\t{}\n".format(getHeader(k, headerKeys), k, bestMiRnasDict[k][p]))
                                out1.write(getStEnDescr(bestMiRnasDict[k][p], headerKeys))
                                out2.write(">" + getHeader(k, headerKeys) + '|' + k + '|' + p + "\n")
                                out2.write(bestMiRnasDict[k][p].split("\t")[0].upper() + "\n")
                                out3.write(">" + getHeader(k, headerKeys) + '|' + k + '|' + p + "\n")
                                out3.write(bestMiRnasDict[k][p].split("\t")[2] + "\n")

def getRange(item):
    s = int(item.split("\t")[0].split("-")[0])
    e = int(item.split("\t")[0].split("-")[1])
    return range(s, e)

def checkIntesects(inList): 
    # check intersects and remove intersected records with lower SVM_score
	delList = []
	for i in range(len(inList)):
		if inList[i] == inList[-1]:
			break
		else:
			aSet = set(getRange(inList[i]))
			for j in range(i+1, len(inList)):
				b = getRange(inList[j])
				if aSet.intersection(b):
					delList.append(inList[j])
	delList = list(set(delList))
	[inList.remove(x) for x in delList]
	return inList

def send2check(fineDict):
    # TE991 4904-4989	3P	uaguaaaacucaggagggaauccgaguaagggauccugcaacuagguucagagguaUUGCUCUCUUCUUUCUCUAAGuuugugcaa	TE991:4904-4989:57_21	UUGCUCUCUUCUUUCUCUAAG	3P	0.903257
    for k in sorted(fineDict.keys()):
        if len(fineDict[k]) > 1:
            fineDict[k].sort(key=lambda x: x.split("\t")[6], reverse=True)  # sort list by SVM_score in descending order
            fineDict[k] = checkIntesects(fineDict[k])
    return fineDict            

def filterBestMiRnaDict(bestMiRnasDict):
    fineDict = {}
    for k in bestMiRnasDict.keys():
        for p in sorted(bestMiRnasDict[k].keys()):
            if bestMiRnasDict[k][p]:
                if k.split(":")[0] not in fineDict.keys():
                    fineDict[k.split(":")[0]] = {}
                    fineDict[k.split(":")[0]] = []
                    fineDict[k.split(":")[0]].append(k.split(":")[1] + "\t" + p + "\t" + bestMiRnasDict[k][p])
                else:
                    fineDict[k.split(":")[0]].append(k.split(":")[1] + "\t" + p + "\t" + bestMiRnasDict[k][p])
    
    fineDictFilt = send2check(fineDict)
    return fineDictFilt

def reformFineDict(fineDictFilt):
    reformFineDict = {}
    for k in sorted(fineDictFilt.keys()):
        for rec in fineDictFilt[k]:
            kOut1 = k + ":" + rec.split("\t")[0]
            kOut2 = rec.split("\t")[1]
            if kOut1 not in reformFineDict.keys():
                reformFineDict[kOut1] = {}
                reformFineDict[kOut1][kOut2] = "\t".join(rec.split("\t")[2:])                
            else:
                if kOut2 not in reformFineDict[kOut1].keys():
                    reformFineDict[kOut1][kOut2] = "\t".join(rec.split("\t")[2:])                    
                else:
                    reformFineDict[kOut1][kOut2] = "\t".join(rec.split("\t")[2:])
    return reformFineDict

def main(miRParaOut, headerKeys):
    inDict = getMirParaDict(miRParaOut) # upload miRPara output
    bestMiRnasDict = filterCandidates(inDict) # filter miRNA with the best SVM_score
    fineDict = filterBestMiRnaDict(bestMiRnasDict)  # here I had to add next filter for pre-miRNAs intersects in order to chose just one with higher SVM_score
    outDict = reformFineDict(fineDict)
    cnt = 0
    for k in outDict.keys():
        for p in outDict[k]:
            if outDict[k][p]:
                cnt += 1
    print("Nr. of filtered miRNAs: {}".format(cnt))
    addHeaderOutput(outDict, headerKeys, miRParaOut) # print results to the text file with original headers in the 1st column: 'Acom|LG01|TE_132|gypsy_Tekay|noTSD'

if __name__ == '__main__':
    main(miRParaOut, headerKeys)
