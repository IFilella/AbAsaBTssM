import glob
import numpy
import fastaf
import matplotlib.pyplot as plt

genestinf = 500
genestsup = 8000
taxids = glob.glob("genomesTssM/*")
count = 0
tssmlist = numpy.genfromtxt("TssM/TssM.list.txt",delimiter=";",skip_header=1,dtype=str)
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
tssjs = tssmlist[:,-1]
tssjs = tssjs.astype(int)
asabs = tssmlist[:,-2]
asabs = asabs.astype(int)
tssms = tssmlist[:,-3]
tssms = tssms.astype(int)
total = len(tssjs)
totalj = sum(tssjs)
totalm = sum(tssms)
totalasab = sum(asabs)
print(total,totalj,totalm,totalasab)
counts = [0,0,0,0]
for line in tssmlist:
    if line[2] == "1" and line[3] == "1": counts[0]+=1
    elif line[2] == "0" and line[3] == "1": counts[1]+=1
    elif line[2] == "1" and line[3] == "0": counts[2]+=1
    elif line[2] == "0" and line[3] == "0": counts[3]+=1
print(counts)
print(sum(counts))
"""
#Get average length of TssM sequences
aux = []
for seq in tssmseqs.seqs:
    print(len(seq))
    aux.append(len(seq))
aux = numpy.asarray(aux)
print(numpy.mean(aux))
