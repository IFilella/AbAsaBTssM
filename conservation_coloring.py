import fastaf
import argparse
import pickle as pickle
import copy

pallette1 = ["#063264","#256BAE","#4594C4","#99C8E0","#CEE4EF","#F2F5F6","#FBD1B9","#F19D7C","#D76350","#B41C2D","#6C021F","#6C021F"]
pallette2 = ["#F9F9FF","#F9F9FF","#D0D0FF","#D0D0FF","#A2A2FF","#A2A2FF","#7372FF","#7372FF","#4545FF","#4545FF","#0D00FF","#0D00FF"]

parser = argparse.ArgumentParser(
    description = '')
parser.add_argument(
    '-cons',
    dest="cons",
    help = "File with conservation values")
parser.add_argument(
    '-seq',
    dest="seq",
    help = "Aligned seq")
parser.add_argument(
    '-mch',
    dest="mch",
    help = "Affected model-chains")
inputfiles = parser.parse_args()
cons=inputfiles.cons
seq=inputfiles.seq
mch=inputfiles.mch
f=open(cons)
cons=[l for l in f]
cons=cons[0].split(",")
cons=cons[1:-1]
cons=[int(float(i)) for i in cons]
f.close()
seq = fastaf.Alignment(seq).seqs[0]
if len(cons)!=len(seq): raise ValueError("The conservation file and the aligned sequence have different length")
rescount=0
for j,res in enumerate(seq):
    if res=="-": continue
    rescount+=1
    #print(rescount,res,cons[j])
    print("color %s,a,r #%s:%d"%(pallette2[cons[j]],mch,rescount))
