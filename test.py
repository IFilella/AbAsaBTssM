import fastaf

#Input data
prot = "TssM"
fafile = fastaf.fastaf("/work/ifilella/AbAsaBTssM/TssM/Acinetobacter.%s.3.fa"%prot)

#Output data
f = open("TssM/NOBaumannii.%s.3.fa"%prot,"w")

count = 0
for homol in fafile.homolseqs:
    title = homol.title
    seq = homol.seq
    if "Baumanni" in title or "baumanni" in title: continue
    else:
        f.write(">" + title+"\n")
        f.write(seq+"\n")
        count+=1
print(count)

"""
fout = open("batchlist.txt","w")
f = open("TssM/TssMctd.complete.aln")
for line in f:
    if line[0]==">" and "TssM" not in line:
        taxid = line.split("|")[-1].replace("\n","")
        print(taxid)
        fout.write(taxid+"\n")
"""
