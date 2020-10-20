#!/usr/bin/env python3

#####
# takze beh progrmau budu tak ze:
# 1. se pujde pres nonNestetTes.fa s puvodnimi headers 'Alyr|scaffold_1|TE_132|copia_Ivana|noTSD'
# 2. checkne jestli ma/nema miRNA (else # 5.)
# 3. if true: u vsech prislusnych miRNA zjisti pozici v TE (LTR/IN) a inkrementuje jejich nr. v ramci TE
# 4. zapise miRNA seq do multifasty s header: '>Alyr|scaffold_1|TE_132|copia_Ivana|noTSD|miRNA:0'
# 5. info zapise do table: spec  chr     teID    fam     tsd miRnaInLTRs miRnaInINs  sum
 
import sys
from Bio import SeqIO

allTesFa = sys.argv[1]
filtMiRParaBed = sys.argv[2]
protGff = sys.argv[3]

def getMiRnaHeadList(filtMiRParaBed):
    miRnaHeadList = []
    with open(filtMiRParaBed) as bed:
        for miRna in bed:
            if miRna.split("\t")[0] not in miRnaHeadList:
                miRnaHeadList.append(miRna.split("\t")[0])
    return miRnaHeadList

def getTesStartsEnds(miRna, protGff):
    headList = miRna.split("\t")[0].split("|")
    spec,chrom,teID,fam,tsd = 0,1,2,3,4
    ch = headList[chrom]
    miRnaStart, miRnaEnd = int(miRna.split("\t")[1]), int(miRna.split("\t")[2])
    teBase = "TE_BASE " + str(headList[teID].split("_")[1])
    teDict = {'miRna_coord':[0,0], 'teCoord':[0,0], 'ltrLeft':[0,0], 'ltrRight':[0,0]}    # initiate teDict
    with open(protGff) as gff:
        for l in gff:                                                                                   # populate teDict
            if ch + "\t" in l and teBase + ";" in l:
                if "\tnested_repeat\t" in l:
                    teDict['teCoord'] = [int(l.split("\t")[3]), int(l.split("\t")[4])]
                    teDict['miRna_coord'] = [teDict['teCoord'][0] + miRnaStart, teDict['teCoord'][0] + miRnaEnd]
                elif "name=ltr left;" in l:
                    teDict['ltrLeft'] = [int(l.split("\t")[3]), int(l.split("\t")[4])]
                elif "name=ltr right;" in l:
                    teDict['ltrRight'] = [int(l.split("\t")[3]), int(l.split("\t")[4])]
    return teDict

def checkIntersects(teMiRNAList, protGff):
    inLtr, inIn = 0, 0
    for miRna in teMiRNAList:  
        teDict = getTesStartsEnds(miRna, protGff)
        if teDict['ltrLeft'][1] < teDict['ltrRight'][0]:
            if teDict['miRna_coord'][0] in range(teDict['ltrLeft'][1], teDict['ltrRight'][0]) and teDict['miRna_coord'][1] in range(teDict['ltrLeft'][1], teDict['ltrRight'][0]):   # IN miRNAs
                inIn += 1
        elif teDict['ltrLeft'][0] > teDict['ltrRight'][1]:
            if teDict['miRna_coord'][0] in range(teDict['ltrRight'][1], teDict['ltrLeft'][0]) and teDict['miRna_coord'][1] in range(teDict['ltrRight'][1], teDict['ltrLeft'][0]):   # IN miRNAs
                inIn += 1
        if teDict['miRna_coord'][0] in range(teDict['ltrLeft'][0], teDict['ltrLeft'][1]) and teDict['miRna_coord'][1] in range(teDict['ltrLeft'][0], teDict['ltrLeft'][1]):   # LTR left miRNAs
            inLtr += 1
        elif teDict['miRna_coord'][0] in range(teDict['ltrRight'][0], teDict['ltrRight'][1]) and teDict['miRna_coord'][1] in range(teDict['ltrRight'][0], teDict['ltrRight'][1]):   # LTR right miRNAs
            inLtr += 1
    return inLtr, inIn

def getLtrInMiRnaIntersects(teHead, filtMiRParaBed, protGff):
    inLtr, inIn = 0, 0
    miRnaHeadList = getMiRnaHeadList(filtMiRParaBed)

    if teHead in miRnaHeadList:
        teMiRNAList = []
        with open(filtMiRParaBed) as bed:
            for miRna in bed:
                if teHead in miRna:
                    teMiRNAList.append(miRna.rstrip())
        return checkIntersects(teMiRNAList, protGff)
    else:
        return inLtr, inIn

def walkThroughFa(allTesFa, filtMiRParaBed, protGff):
    outDict = {}
    for te in SeqIO.parse(allTesFa, "fasta"):
        inLtr, inIn = getLtrInMiRnaIntersects(te.id, filtMiRParaBed, protGff)
        outDict[te.id] = {'LTR': inLtr, 'IN': inIn}

    outTabName = allTesFa.replace(".fa","_miRnas.tsv")
    with open(outTabName, "w") as out:
        for k in sorted(outDict.keys()):
            out.write("{}\t{}\t{}\t{}\n".format("\t".join(k.split("|")), str(outDict[k]['LTR']), str(outDict[k]['IN']), str(outDict[k]['LTR'] + outDict[k]['IN'])))

def main(allTesFa, filtMiRParaBed, protGff):
    walkThroughFa(allTesFa, filtMiRParaBed, protGff)

if __name__ == "__main__":
    main(allTesFa, filtMiRParaBed, protGff)