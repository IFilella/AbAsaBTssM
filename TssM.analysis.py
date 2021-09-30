import glob
import numpy
import fastaf
import matplotlib.pyplot as plt
from matplotlib.dates import date2num
import datetime
import closeness
import numpy as np

#Input Data
genestinf = 500
genestsup = 8000
taxids = glob.glob("genomesTssM/*")
count = 0
tssmlist = numpy.genfromtxt("TssM/TssM.list.2.txt",delimiter=";",skip_header=1,dtype=str)
tssmseqs = fastaf.Alignment("/work/ifilella/AbAsaBTssM/TssM/Acinetobacter.TssM.fasta")
"""
#Get total genes per genome
totalgenes = []
totaltaxs = []
for taxid in taxids:
    tax = taxid.split("/")[1]
    taxgenome = fastaf.fastaf("%s/%s.fasta"%(taxid,tax))
    taxgenes = len(taxgenome.homolseqs)
    if taxgenes > genestinf and taxgenes < genestsup:
        count+=1
        totalgenes.append(taxgenes)
        totaltaxs.append(tax)

totalgenes, totaltaxs = zip(*sorted(zip(totalgenes, totaltaxs)))
print(list(zip(totalgenes,totaltaxs)))
print(len(totalgenes))
print(numpy.mean(totalgenes))
print(numpy.std(totalgenes))

#Plot total genes per genome
totalgenes = numpy.asarray(totalgenes)
plt.figure()
plt.hist(totalgenes, bins=25,density=False,rwidth=0.93)
plt.xlabel("Number of genes")
plt.ylabel("Frequency")
plt.savefig("Acinetobacter_hist.png")
plt.show()

#Get % of TssM,TssJ and AsaB
tssms = tssmlist[:,2]
tssms = tssms.astype(int)
asabs = tssmlist[:,3]
asabs = asabs.astype(int)
tssjs = tssmlist[:,4]
tssjs = tssjs.astype(int)
tssbs = tssmlist[:,5]
tssbs = tssbs.astype(int)
tssks = tssmlist[:,6]
tssks = tssks.astype(int)
total = len(tssms)
totalm = sum(tssms)
totalasab = sum(asabs)
totalj = sum(tssjs)
totalb = sum(tssbs)
totalk = sum(tssks)
print(total,totalm,totalasab,totalj,totalb,totalk)
counts = [0,0,0,0]
for line in tssmlist:
    if line[5] == "0": continue
    if line[2] == "1" and line[3] == "1": counts[0]+=1
    elif line[2] == "0" and line[3] == "1": counts[1]+=1
    elif line[2] == "1" and line[3] == "0": counts[2]+=1
    elif line[2] == "0" and line[3] == "0": counts[3]+=1
print(counts)
print(sum(counts))

#Get average length of TssM sequences
aux = []
for seq in tssmseqs.seqs:
    print(len(seq))
    aux.append(len(seq))
aux = numpy.asarray(aux)
print(numpy.mean(aux))

# Get the conservation level of M, AsaB, B and K for A.b and Acinetonacter spp
aliAciM = fastaf.Alignment(aliname="/work/ifilella/AbAsaBTssM/TssM/Acinetobacter.TssM.2.aln")
aliAciA = fastaf.Alignment(aliname="/work/ifilella/AbAsaBTssM/TssM/Acinetobacter.AsaB.2.aln")
aliAciB = fastaf.Alignment(aliname="/work/ifilella/AbAsaBTssM/TssM/Acinetobacter.TssB.2.aln")
aliAciK = fastaf.Alignment(aliname="/work/ifilella/AbAsaBTssM/TssM/Acinetobacter.TssK.2.aln")
aliAbM = fastaf.Alignment(aliname="/work/ifilella/AbAsaBTssM/TssM/Baumannii.TssM.2.aln")
aliAbA = fastaf.Alignment(aliname="/work/ifilella/AbAsaBTssM/TssM/Baumannii.AsaB.2.aln")
aliAbB = fastaf.Alignment(aliname="/work/ifilella/AbAsaBTssM/TssM/Baumannii.TssB.2.aln")
aliAbK = fastaf.Alignment(aliname="/work/ifilella/AbAsaBTssM/TssM/Baumannii.TssK.2.aln")

x = ["TssM","AsaB","TssB","TssK"]
y = []
z = []
legend = ["Acinetonacter","A. Baumannii"]
N = len(x)
ind = np.arange(N)
width = 0.37

for i,ali in enumerate([aliAciM,aliAbM,aliAciA,aliAbA,aliAciB,aliAbB,aliAciK,aliAbK]):
    seqs = ali.seqs
    seqsvec = closeness._seqs2vec(seqs)
    random_score = closeness.get_random_score(seqs,seqs,nsample=10)
    identical_score = closeness.get_identical_score(seqs,seqs)
    score = closeness.get_subscore_mixvec(seqsvec,seqsvec)
    print(random_score,identical_score,score)
    similarity = ((float(score)-float(random_score))/(float(identical_score)-float(random_score)))
    print(similarity)
    if i%2 == 0: y.append(similarity)
    else: z.append(similarity)
print(x,y,z)

fig = plt.figure()
ax = fig.add_subplot(111)
rects1 = ax.bar(ind, y, width)
rects2 = ax.bar(ind+width, z, width)
ax.set_ylabel('Conservation Level')
ax.set_xticks(ind+width/len(legend))
ax.set_xticklabels(x)
ax.legend( (rects1[0], rects2[0]), legend, bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=1, ncol=2)

def autolabel(rects):
    for rect in rects:
        h = rect.get_height()
        ax.text(rect.get_x()+rect.get_width()/2., 1.01*h, '%.3f'%float(h),
                ha='center', va='bottom')

autolabel(rects1)
autolabel(rects2)
plt.savefig("conservation.png")
plt.savefig("conservation.pdf")
plt.show()
"""
