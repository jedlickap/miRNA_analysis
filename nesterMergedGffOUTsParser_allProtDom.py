#!/usr/bin/env python3

import sys
import re
from operator import itemgetter
from Bio import SeqIO
import subprocess
import os

headerKeys = sys.argv[1]
gff = sys.argv[2]
print(headerKeys)
with open(headerKeys) as headKeys:
    for rec in headKeys:
        print(rec)
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

def getRtGff(gff, inDict, headerKeys):
    cnt = 0
    with open(headerKeys) as headKeys:
        for rec in headKeys:
            recTrue = False
            ch, teID = rec.split("\t")[0].split("|")[1], rec.split("\t")[0].split("|")[2].split("_")[1]
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
                    print(inDict[ch][teID])
                    print(domDict)
                    recTrue = True
                domDict['name=CHD'] = True
                if False not in list(domDict.values()):
                    print(inDict[ch][teID])
                    print(domDict)
                    recTrue = True
                if recTrue:
                    cnt +=1
    print("nr of TEs with all protDoms: {}".format(cnt))

    #     with open(gff) as file:
    #         for l in file:
    #             if "##" in l:
    #                 out.write(l)
    #     for ch in sorted(inDict.keys()):
    #         l = len(inDict[ch].keys())
    #         for i in range(0, l):
    #             #if "polypeptide_conserved_region" in inDict[ch][str(i)].keys() and "target_site_duplication" in inDict[ch][str(i)].keys():  # filter for presence of PCR and TSD
    #             if "polypeptide_conserved_region" in inDict[ch][str(i)].keys():     # filter for presence of PCR only
    #                 rt = False
    #                 for dom in inDict[ch][str(i)]["polypeptide_conserved_region"]:
    #                     if "name=RT" in dom[2]:
    #                         annotRT = dom[2].split(";")[2].split(" ")[0]
    #                         rt =True
    #                     else:
    #                         annotElse = dom[2].split(";")[2].split(" ")[0] 
    #                 if rt:
    #                     annot = annotRT
    #                 else:
    #                     annot = annotElse
    #                 with open(gff) as file:
    #                     for l in file:
    #                         if ch in l and "nested_repeat" in l and "TE_BASE " + str(i) + ";" in l:
    #                             lstart = "\t".join(l.rstrip().split("\t")[0:-1])
    #                             lend = l.rstrip().split("\t")[8].split(";")[0] + ";" + l.rstrip().split("\t")[8].split(";")[1] + ";" + annot + ";" + l.rstrip().split("\t")[8].split(";")[-1]
    #                             l = lstart + "\t" + lend
    #                             out.write(l + "\n")
    #                         if ch in l and "TE_BASE " + str(i) + ";" in l and "nested_repeat" not in l:
    #                             out.write(l)
    # return gff.replace(".gff", "_protDom.gff")

inDict = fillInDict(gff)
getRtGff(gff, inDict, headerKeys)

def getAllTesGff(gff, inDict):
    annot = ""
    with open(gff.replace(".gff", "_allTEs.gff"), "w") as out:     # alternativelly "RT" domain
        with open(gff) as file:
            for l in file:
                if "##" in l:
                    out.write(l)
        for ch in sorted(inDict.keys()):
            l = len(inDict[ch].keys())
            for i in range(0, l):
                with open(gff) as file:
                    for l in file:
                        if ch in l and "nested_repeat" in l and "TE_BASE " + str(i) + ";" in l:
                            lstart = "\t".join(l.rstrip().split("\t")[0:-1])
                            lend = l.rstrip().split("\t")[8].split(";")[0] + ";" + l.rstrip().split("\t")[8].split(";")[1] + ";" + annot + ";" + l.rstrip().split("\t")[8].split(";")[-1]
                            l = lstart + "\t" + lend
                            out.write(l + "\n")
                        if ch in l and "TE_BASE " + str(i) + ";" in l and "nested_repeat" not in l:
                            out.write(l)
    return gff.replace(".gff", "_allTEs.gff")

def getTeStartEnd(recList):
    return int(recList[3]), int(recList[4])

def getRfCoordList(sortListRec):
    rfCoordList = []
    for rfRec in sortListRec[6]:
        rfCoordList.append(rfRec[0])
        rfCoordList.append(rfRec[1])
    return rfCoordList

def getNestedNonnestedTEs(inDict):              # list structure is: [teLen, chrom, teID, start, end, [repFragList]] !!!for one chromosome!!!
    N_NnTesList = []
    for ch in sorted(inDict.keys()):
        for te in inDict[ch].keys():
            teS = inDict[ch][te]["nested_repeat"][0]
            teE = inDict[ch][te]["nested_repeat"][1]
            teLen = int(inDict[ch][te]["nested_repeat"][1]) - int(inDict[ch][te]["nested_repeat"][0])
            rfList = []
            if "annot" in inDict[ch][te]["nested_repeat"][2].split(";")[2]:
                annot = inDict[ch][te]["nested_repeat"][2].split(";")[2].replace("annot=", "")
            else:
                annot = "no-annotation"
            for rf in inDict[ch][te]["repeat_fragment"]:
                    rfList.append(rf)                
            toAppend = [str(teLen), ch, "TE_BASE " + te, teS, teE, annot, rfList]
            N_NnTesList.append(toAppend)
    sortList =sorted(N_NnTesList, key = lambda x: int(x[0]))
    allList = []
    nestList = []
    nonNestList = []
    allOrigList = []
    teCnt = len(sortList)
    for i in range(teCnt):
        shortTE_s, shortTE_e = getTeStartEnd(sortList[i])
        origList = []            
        if i < teCnt:
            for j in range(i + 1, teCnt):
                rfCoordList = getRfCoordList(sortList[j])
                longTE_s, longTE_e = getTeStartEnd(sortList[j])
                if shortTE_s > longTE_s and shortTE_e < longTE_e:
                    allOrigList.append("\t".join(sortList[j][1:-1]))
                    if str(shortTE_s) in rfCoordList and str(shortTE_e) in rfCoordList:     # here, the coordinates of nested TE are checked with repeat fragments coordinates of original TE
                        origList.append(sortList[j][1:-1])                        
        toAppend = [sortList[i][1:-1], "IN", origList]
        allList.append(toAppend)
    
    allOrigList_uniq = list(set(allOrigList))
        
    for i in range(len(allList)):
        if allList[i][2]:
            nestList.append(allList[i])
        else:
            ch = allList[i][0][0]
            teID =allList[i][0][1]
            nonNest = True
            for rec in allOrigList_uniq:
                if ch + "\t" + teID + "\t" in rec:
                    nonNest =False
            if nonNest:
                nonNestList.append(allList[i])
    return nonNestList, nestList

def getltrList(ch, ID, gff):
    ltrLeftList = []
    with open(gff) as f:
        for l in f:
            if ch + "\t" in l and ID + ";" in l and "ID=LTR LEFT" in l:
                ltrLeftList.append(l.rstrip())
    ltrRightList = []
    with open(gff) as f:
        for l in f:
            if ch + "\t" in l and ID + ";" in l and "ID=LTR RIGH" in l:
                ltrRightList.append(l.rstrip())
    return [ltrLeftList, ltrRightList]

def getLtrFaIdent(ltrLists, genomeFa):
    # getBed
    ltrLenList = [0, 0]
    filesList = []
    for i in range(len(ltrLists)):
        if len(ltrLists[i]) == 1:
            for r in ltrLists[i]:
                r = r.split("\t")
                with open("ltr.bed", "w") as out1:
                    out1.write(r[0] + "\t" + r[3] + "\t" + r[4] + "\n")
                    ltrLen = int(r[4]) - int(r[3])
                command = 'bedtools getfasta -fi {} -bed {} -fo ltr{}.fa'.format(genomeFa, "ltr.bed", i)
                process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
                process.wait()
                print(process.returncode)
                ltrLenList[i] = ltrLen
                filesList.append("ltr.bed")
                filesList.append("ltr{}.fa".format(i))
        elif len(ltrLists[i]) > 1:
            with open("ltr.bed", "w") as out1:                    
                nrFrag = len(ltrLists[i])
                for j in range(nrFrag):
                    r = ltrLists[i][j].split("\t")
                    out1.write(r[0] + "\t" + r[3] + "\t" + r[4] + "\n")
                    
            faName = "ltr{}_fragm.fa".format(i)
            command = 'bedtools getfasta -fi {} -bed {} -fo {}'.format(genomeFa, "ltr.bed", faName)
            process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
            process.wait()
            print(process.returncode)
            seq = ""
            with open("ltr{}.fa".format(i), "w") as faOut:
                faOut.write(">ltr_joined_{}\n".format(i))
                for rec in SeqIO.parse(faName, "fasta"):
                    seq += str(rec.seq)
                faOut.write(str(seq) + "\n")
            ltrLen = len(seq)    
            ltrLenList[i] = ltrLen
            filesList.append("ltr.bed")
            filesList.append(faName)
    # run stretcher
    if ltrLenList[0] > 0 and ltrLenList[1] > 0:
        cmd = "/home/pavel/Documents/Soft/EMBOSS-6.6.0/emboss/stretcher " + "ltr0.fa ltr1.fa -outfile ltrIdent.txt"
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        process.wait()
        print(process.returncode)        
    identReg = r' \((\d+\.\d+)%'
    ident = 0
    avgLtrLen = (int(ltrLenList[0]) + int(ltrLenList[1])) / 2
    with open("ltrIdent.txt") as identFile:
        for l in identFile:
            if "# Identity:" in l:
                ident = float(re.search(identReg, l.rstrip()).group(1))
    filesList.append("ltrIdent.txt")

    # get Kimura distance:
    ## run clustalw to get .phy file
    k80 = ""
    if not os.path.exists("outfile"):
        with open("outfile", "w") as of:
            of.write("")
    if ltrLenList[0] > 0 and ltrLenList[1] > 0:
        # merge two fa files:
        cmd = "cat ltr0.fa ltr1.fa >> infile.fa"
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        process.wait()
        
        cmd1 = "clustalw infile.fa -output=PHYLIP"
        process = subprocess.Popen(cmd1, shell=True, stdout=subprocess.PIPE)
        process.wait()
        print(process.returncode)

        os.rename('infile.phy', 'infile')
        with open("opt_dnadist", "w") as out:
            out.write("infile\nR\nD\nY\n")

        cmd1 = "cat opt_dnadist | phylip dnadist"
        process = subprocess.Popen(cmd1, shell=True, stdout=subprocess.PIPE)
        process.wait()
        print(process.returncode)

        with open("outfile") as outfile:
            next(outfile)
            for l in outfile:
                lList = l.rstrip().split(" ")
                k80 = lList[2]  # return K80
                break
        filesList.append("infile")
        filesList.append("infile.fa")
        filesList.append("infile.dnd")
        
    if not k80:
        k80 = "NA"        
    
    # remove redundant files
    for file in filesList:
        if os.path.exists(file):
            os.remove(file)                
    return [str(ident), str(avgLtrLen), str(ltrLenList[0]), str(ltrLenList[1]), str(k80)]           

def getNestOrigLTRident(nestList, nonNestList, gff, genomeFa):
    outNestOrigList = []
    outAllTesList = []
    for pair in nestList:
        nestList = []
        origList = []
        ch = pair[0][0]
        nestAnnot, origAnnot = "", ""
        nestList.append(pair[0][1])
        if "::" in pair[0][-1]:
            nestAnnot = pair[0][-1].split("::")[2].split("/")[1] + "_" + pair[0][-1].split("::")[-1]
        else:
            nestAnnot = "no-annotation"
        nestList.append(nestAnnot)
        
        origList.append(pair[2][0][1])
        if "::" in pair[2][0][-1]:
            origAnnot = pair[2][0][-1].split("::")[2].split("/")[1] + "_" + pair[2][0][-1].split("::")[-1]
        else:
            origAnnot = "no-annotation"
        origList.append(origAnnot)
        
        nestList.append("\t".join(getLtrFaIdent(getltrList(ch, nestList[0], gff), genomeFa)))
        origList.append("\t".join(getLtrFaIdent(getltrList(ch, origList[0], gff), genomeFa)))
        
        nestStr = "\t".join(nestList)
        origStr = "\t".join(origList)
        outNestOrigList.append(ch + "\t" + nestStr + "\t" + origStr)
        outAllTesList.append(ch + "\t" + nestStr + "\tnested")
        outAllTesList.append(ch + "\t" + origStr + "\toriginal")
    with open("{}_nestOrigIdent.tsv".format(gff.replace(".gff", "")), "w") as out:
        for l in outNestOrigList:
            out.write(l + "\n")
    for te in nonNestList:
        ch = te[0][0]
        ID = te[0][1]
        if "::" in te[0][4]:
            annot = te[0][4].split("::")[2].split("/")[1] + "_" + te[0][4].split("::")[-1]
        else:
            annot = "no-annotation"
        outAllTesList.append(ch + "\t" + ID + "\t" + annot + "\t" + "\t".join(getLtrFaIdent(getltrList(ch, ID, gff), genomeFa)) + "\tnon-nested")
    with open("{}_allTeLtrIdent.tsv".format(gff.replace(".gff", "")), "w") as out:
        for l in outAllTesList:
            out.write(l + "\n")
        


# for k in inDict.keys():
#     for k1 in inDict[k]:
#         for k2 in inDict[k][k1]:
#             print(k, k1, k2, inDict[k][k1][k2])

#nonNestList1, nestList1 = getNestedNonnestedTEs(inDict)
#getNestOrigLTRident(nestList1, nonNestList1, gff, genomeFa)

# rtGffName = getRtGff(gff, inDict)
# rtInDict = fillInDict(rtGffName)
# nonNestList2, nestList2 = getNestedNonnestedTEs(rtInDict)
#getNestOrigLTRident(nestList2, nonNestList2, rtGffName, genomeFa)
