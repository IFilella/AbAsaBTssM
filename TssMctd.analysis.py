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
gxfile = numpy.genfromtxt("TssMctd/GXXXGXXXG.noAcineto.2.txt",dtype=str,skip_header=1,delimiter=";")
total = len(gxfile)

#Output data
totalb = 0
totalp = 0
totalnp = 0
totalp1 = 0
totalp2 = 0
totalp3 = 0
totalplus = 0
Gtssj = 0
NOGtssj = 0

for line in gxfile:
    line = line.tolist()
    if "Error" in line[-1]: continue
    if line[-2] == "0": continue
    totalb+=1
    if line[0] != "NO":
        totalp +=1
        patterns = line[0].split("_")
        if len(patterns)==1: totalp1+=1
        if len(patterns)==2: totalp2+=1
        if len(patterns)==3: totalp3+=1
        if len(patterns)>3: totalplus+=1
        if line[-1] == "1": Gtssj+=1
    else:
        totalnp+=1
        if line[-1] == "1":
            NOGtssj+=1

print("Total: %s\nFilterTssB: %s\nGxxxGxxxG: %s\np1: %s\np2: %s\np3: %s\np+: %s\nNOGxxxGxxxG: %s"%(total,totalb,totalp,totalp1,totalp2,totalp3,totalplus,totalnp))
print("GxxxGxxxG with TssJ: %s\nNO GxxxGxxxG with TssJ: %s"%(Gtssj,NOGtssj))


###Find GXXGXXXG pattern on TssMctd of Acinetobacter species###

#Input data
Blist = numpy.genfromtxt("TssM/TssM.list.3.txt",skip_header=1,delimiter=";",dtype=str)
alictdAci = fastaf.Alignment(aliname="TssM/Acinetobacter.TssM.3.fa")
alictdBau = fastaf.Alignment(aliname="TssM/Baumannii.TssM.3.fa")
alictdNOBau = fastaf.Alignment(aliname="TssM/NOBaumannii.TssM.3.fa")
alisctd = [alictdAci,alictdBau,alictdNOBau]
alisnames = ["Acinetobacter","Baumannii","NOBaumannii"]

#Output data
fA = open("TssMctd/GXXXGXXXG.Acineto.3.txt","w")
fB = open("TssMctd/GXXXGXXXG.Baumannii.3.txt","w")
fN = open("TssMctd/GXXXGXXXG.NOBaumannii.3.txt","w")
files = [fA,fB,fN]

for j,alictd in enumerate(alisctd):
    print("---------------------------------------------")
    print(j,alisnames[j])
    total = 0
    totalb = 0
    totalp = 0
    totalp1 = 0
    totalp2 = 0
    totalp3 = 0
    totalplus = 0
    totalnp = 0
    for i,name in enumerate(alictd.names):
        total+=1
        unaliseq = alictd.seqs[i]
        taxid = name.split("_")[-1]
        organism = "_".join(name.split("_")[0:-2])
        indx = numpy.where(Blist==taxid)[0][0]
        presenceB = list(Blist[indx])[5]
        if presenceB != '1': continue
        totalb+=1
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
        else:
            totalnp+=1
            files[j].write("NOTFOUND;%s;%s\n"%(organism,taxid))
    print("Total: %s\nFilterTssB: %s\nGxxxGxxxG: %s\np1: %s\np2: %s\np3: %s\np+: %s\nNOGxxxGxxxG: %s"%(total,totalb,totalp,totalp1,totalp2,totalp3,totalplus,totalnp))

