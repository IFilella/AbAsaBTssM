import fastaf
import argparse
import pickle as pickle
import copy

pallette1 = ["#F9F9FF","#F9F9FF","#D0D0FF","#D0D0FF","#A2A2FF","#A2A2FF","#7372FF","#7372FF","#4545FF","#4545FF","#0D00FF","#0D00FF"]
pallette2 = ["#FF0000","#FFA300","#FFFF00","#00FF00","#00FFFF","#0000FF"]
pallette3 = ["#FF0000","#FF0000","#FFA300","#FFA300","#FFFF00","#FFFF00","#00FF00","#00FF00","#00FFFF","#00FFFF","#0000FF","#0000FF"]
pallette4 = ["#FF0000","#FF0000","#FF0000","#FF0000","#00FF00","#00FF00","#00FF00","#00FFFF","#0000FF","#0000FF","#0000FF","#0000FF"]

parser = argparse.ArgumentParser(
    description = '')
parser.add_argument(
    '-data',
    dest="data",
    help = "File with residue conservation (JalView) or residue lDDT confidence (AlphaFold2)")
parser.add_argument(
    '-dtype',
    dest = "dtype",
    help = "Data type, either 'cons' for conservation  or 'conf' for lDDT confidence")
parser.add_argument(
    '-seq',
    dest="seq",
    help = "Aligned seq when conservation or regular seq when confidence")
parser.add_argument(
    '-mch',
    dest="mch",
    help = "Affected model-chains")

#Parsing input data
inputfiles = parser.parse_args()
data = inputfiles.data
dtype = inputfiles.dtype
if dtype != 'cons' and dtype != 'conf': raise ValueError("dtype must be either 'cons' or 'conf'")
seq=inputfiles.seq
mch=inputfiles.mch

seq = fastaf.Alignment(seq).seqs[0]
print(seq,len(seq))
f = open(data)
data = [l for l in f]
data = data[0].split(",")
if dtype == 'cons':
    con = [int(float(i)) for i in data]
    print(con,len(con))
elif dtype == 'conf':
    con = [float(i) for i in data]
    print(con,len(con))
f.close()
if len(con)!=len(seq): raise ValueError("The conservation file and the aligned sequence have different length")

#Generate the output for Chimera
rescount=0
if dtype == 'cons':
    for j,res in enumerate(seq):
        if res=="-": continue
        rescount+=1
        #print(rescount,res,con[j])
        print("color %s,a,r #%s:%d"%(pallette4[con[j]],mch,rescount))
elif dtype == 'conf':
    for j,res in enumerate(seq):
        rescount+=1
        #print(rescount,res,con[j])
        if con[j] <= 50: c = 0
        if con[j] > 50 and con[j] <= 60: c = 1
        if con[j] > 60 and con[j] <= 70: c = 2
        if con[j] > 70 and con[j] <= 80: c = 3
        if con[j] > 80 and con[j] <= 90: c = 4
        if con[j] > 90: c = 5
        print("color %s,a,r #%s:%d"%(pallette2[c],mch,rescount))
