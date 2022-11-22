import numpy
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib
import matplotlib.pyplot
import matplotlib.ticker
import pylab

#Input data
homologs = numpy.genfromtxt("AsaB/AsaB.txt",delimiter=";",skip_header=1,dtype=str)

#Output data
count=0
cleanhomols = []
evalues = []
coverages = []
labels = []
tssms = []
tssjs = []
identifiers = []

for homol in homologs:
    if "Error" not in homol[-1]:
        cleanhomols.append(homol)
        count+=1
        identifiers.append(homol[0].split("_")[1])
        evalues.append(float(homol[1]))
        coverages.append(float(homol[2]))
        tssms.append(int(homol[5]))
        tssjs.append(int(homol[6]))
        if "Acineto" in homol[4]:
            labels.append("Acinetobacter")
        elif "Pseudo" in homol[4]:
            labels.append("Pseudomonas")
        elif "Xantho" in homol[4]:
            labels.append("Xanthomonas")
        elif "Achromo" in homol[4]:
            labels.append("Achromobacter")
        else:
            labels.append("Others")

print(len(homologs),count)
print(len(evalues),len(coverages),len(labels))
print(tssms)
print(tssjs)
tssms = list(map(bool,tssms))
tssjs = list(map(bool,tssjs))

dfhomols = pd.DataFrame({"Identifier": identifiers, "Evalue": evalues, "Coverage": coverages, "Phylum": labels, "TssM": tssms, "TssJ": tssjs})
dfhomolsg = dfhomols.groupby("Phylum")
print("Total number of AsaB with a unique genome associated: " + str(len(dfhomols)))
dfacineto = dfhomols.loc[dfhomols['Phylum'] == "Acinetobacter"]
dfnoacineto = dfhomols.loc[dfhomols['Phylum'] != "Acinetobacter"]
print("Acinetobacter genomes with an AsaB homolog: " + str(len(dfacineto)))
print("Acinetobacter genomes with an AsaB sequence: " + str(len(dfacineto.loc[(dfacineto["Evalue"] <= 1e-20) & (dfacineto["Coverage"] >= 0.6)])))
print("Acinetobacter genomes with an AsaB sequence and a TssM homolog: " + str(len(dfacineto.loc[(dfacineto["Evalue"] <= 1e-20) & (dfacineto["Coverage"] >= 0.6) & (dfhomols['TssM'] == True)])))
print("Acinetobacter genomes with an AsaB sequence and a TssJ homolog: " + str(len(dfhomols.loc[(dfacineto["Evalue"] <= 1e-20) & (dfacineto["Coverage"] >= 0.6) & (dfhomols['TssJ'] == True)])))
print("Non Acinetobacter genomes with an AsaB homolog: " + str(len(dfnoacineto)))
print("Non Acinetobacter genomes with an AsaB similar: " + str(len(dfnoacineto.loc[(dfnoacineto["Evalue"] > 1e-20) | (dfnoacineto["Coverage"] < 0.6)])))
print("Non Acinetobacter genomes with an AsaB similar a TssM homolog and a TssJ homolog: " + str(len(dfnoacineto.loc[(dfnoacineto['TssM'] == True) & (dfnoacineto['TssJ'] == True) ])))
print("Non Acinetobacter genomes with an AsaB similar a TssM homolog and without TssJ homolog: " + str(len(dfnoacineto.loc[(dfnoacineto['TssM'] == True) & (dfnoacineto['TssJ'] == False) ])))
print("Non Acinetobacter genomes with an AsaB similar without TssM homolog and with a TssJ homolog: " + str(len(dfnoacineto.loc[(dfnoacineto['TssM'] == False) & (dfnoacineto['TssJ'] == True) ])))
print("Non Acinetobacter genomes with an AsaB similar without a TssM homolog and without a TssJ homolog: " + str(len(dfnoacineto.loc[(dfnoacineto['TssM'] == False) & (dfnoacineto['TssJ'] == False) ])))
exit()
plt.figure()
ax = plt.gca()
for name, homol in dfhomolsg:
    ax.plot(homol["Evalue"], homol["Coverage"], marker="o", linestyle="", label=name,alpha=0.4)
ax.set_xscale('log')
ax.set_xlabel('Blastp e-value (logscale)')
ax.set_ylabel('Coverage')
plt.title("AsaB homologs: E-value against Coverage")
plt.legend()
pylab.tight_layout()
plt.savefig("pictures/evaluecoverage.png")
plt.savefig("pictures/evaluecoverage.pdf")
plt.show()
