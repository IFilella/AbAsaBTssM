import fastaf
import numpy
import pandas as pd
import re
import os
import glob
import search
import subprocess

###Analysis of the GxxxGxxxG output for non Acinetobacter###

#Input data
gxfile = numpy.genfromtxt("TssMctd/GXXXGXXXG.noAcineto.txt",dtype=str,skip_header=1,delimiter=";")
total = len(gxfile)

#Output data
totalB = 0
totalp = 0
totalnp = 0
totalp1 = 0
totalp2 = 0
totalp3 = 0
totalplus = 0
totalGA = 0
totalNGA = 0
totalJ = [0,0,0,0,0,0,0,0]

for line in gxfile:
    line = line.tolist()
    if "Error" in line[-1]: continue
    if line[4] == "0": continue
    totalB+=1
    if line[0] != "NO":
        totalp +=1
        patterns = line[0].split("_")
        if len(patterns)==1: totalp1+=1
        if len(patterns)==2: totalp2+=1
        if len(patterns)==3: totalp3+=1
        if len(patterns)>3: totalplus+=1
        if line[6] == "1":
            totalGA+=1
            print(line)
            if line[5] == "1": totalJ[0]+=1
            elif line[5] == "0": totalJ[1]+=1
        elif line[6] == "0":
            if line[5] == "1": totalJ[2]+=1
            elif line[5] == "0": totalJ[3]+=1
    else:
        totalnp+=1
        if line[6] == "1":
            totalNGA+=1
            print(line)
            if line[5] == "1": totalJ[4]+=1
            elif line[5] == "0": totalJ[5]+=1
        elif line[6] == "0":
            if line[5] == "1": totalJ[6]+=1
            elif line[5] == "0": totalJ[7]+=1

print("Total: %s\nFilterTssB: %s\nGxxxGxxxG: %s\np1: %s\np2: %s\np3: %s\np+: %s\nNOGxxxGxxxG: %s"%(total,totalB,totalp,totalp1,totalp2,totalp3,totalplus,totalnp))
print("GxxxGxxxG with AsaB: %s\nNO GxxxGxxxG with AsaB: %s"%(totalGA,totalNGA))
print(totalJ)
exit()
###Find GXXGXXXG pattern on TssMctd of Acinetobacter species###

#Input data
proteinlist = numpy.genfromtxt("TssM/TssM.list.txt",skip_header=1,delimiter=";",dtype=str)
alictdAci = fastaf.Alignment(aliname="TssM/Acinetobacter.TssM.fa")
alictdBau = fastaf.Alignment(aliname="TssM/Baumannii.TssM.fa")
alictdNOBau = fastaf.Alignment(aliname="TssM/NOBaumannii.TssM.fa")
alisctd = [alictdAci,alictdBau,alictdNOBau]
alisnames = ["Acinetobacter","Baumannii","NOBaumannii"]

#Output data
fA = open("TssMctd/GXXXGXXXG.Acineto.txt","w")
fB = open("TssMctd/GXXXGXXXG.Baumannii.txt","w")
fN = open("TssMctd/GXXXGXXXG.NOBaumannii.txt","w")
files = [fA,fB,fN]

for j,alictd in enumerate(alisctd):
    print("---------------------------------------------")
    print(j,alisnames[j])
    total = 0
    totalB = 0
    totalp = 0
    totalp1 = 0
    totalp2 = 0
    totalp3 = 0
    totalplus = 0
    totalnp = 0
    totalGA = 0
    totalNGA = 0
    for i,name in enumerate(alictd.names):
        total+=1
        unaliseq = alictd.seqs[i]
        taxid = name.split("_")[-1]
        organism = "_".join(name.split("_")[0:-2])
        indx = numpy.where(proteinlist==taxid)[0][0]
        presenceB = list(proteinlist[indx])[5]
        presenceA = list(proteinlist[indx])[3]
        if presenceB != '1': continue
        totalB+=1
        #Found the pattern on the sequence
        patterns1 = re.findall('G..G..G',unaliseq[-100:])
        patterns2 = re.findall('G...G...G',unaliseq[-100:])
        patterns3 = re.findall('G..G...G',unaliseq[-100:])
        patterns4 = re.findall('G...G..G',unaliseq[-100:])
        patterns = patterns1 + patterns2 + patterns3 + patterns4
        if len(patterns)>0:
            totalp+=1
            #print(i,name,patterns)
            if len(patterns)==1: totalp1+=1
            if len(patterns)==2: totalp2+=1
            if len(patterns)==3: totalp3+=1
            if len(patterns)>3: totalplus+=1
            #print("%s;%s;%s\n"%('_'.join(patterns),organism,taxid))
            files[j].write("%s;%s;%s\n"%('_'.join(patterns),organism,taxid))
            if presenceA == '1': totalGA+=1
        else:
            totalnp+=1
            if presenceA == '1': totalNGA+=1
            files[j].write("NOTFOUND;%s;%s\n"%(organism,taxid))
    print("Total: %s\nFilterTssB: %s\nGxxxGxxxG: %s\np1: %s\np2: %s\np3: %s\np+: %s\nNOGxxxGxxxG: %s"%(total,totalB,totalp,totalp1,totalp2,totalp3,totalplus,totalnp))
    print("GxxxGxxxG with AsaB: %s\nNO GxxxGxxxG with AsaB: %s"%(totalGA,totalNGA))
