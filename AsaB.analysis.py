import numpy
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib
import matplotlib.pyplot
import matplotlib.ticker
import pylab

homologs = numpy.genfromtxt("AsaB/AsaB.txt",delimiter=";",skip_header=1,dtype=str)

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
print(dfhomols)
print(dfhomols.loc[(dfhomols['Phylum'] == "Acinetobacter") ])
print(dfhomols.loc[(dfhomols['Phylum'] == "Acinetobacter") & (dfhomols['TssM'] == True)])
print(dfhomols.loc[(dfhomols['Phylum'] == "Acinetobacter") & (dfhomols['TssJ'] == True)])
print(dfhomols.loc[(dfhomols['Phylum'] != "Acinetobacter") ])
print(dfhomols.loc[(dfhomols['Phylum'] != "Acinetobacter") & (dfhomols['TssM'] == True)])
print(dfhomols.loc[(dfhomols['Phylum'] != "Acinetobacter") & (dfhomols['TssJ'] == True)])
plt.figure()
ax = plt.gca()
for name, homol in dfhomolsg:
    ax.plot(homol["Evalue"], homol["Coverage"], marker="o", linestyle="", label=name,alpha=0.4)
ax.set_xscale('log')
#x_major = matplotlib.ticker.LogLocator(base = 10.0, numticks = 5)
#ax.xaxis.set_major_locator(x_major)
#x_minor = matplotlib.ticker.LogLocator(base = 10.0, numticks = 10)
#ax.xaxis.set_minor_locator(x_minor)
#ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_xlabel('Blastp e-value (logscale)')
ax.set_ylabel('Coverage')
plt.title("AsaB homologs: E-value against Coverage")
plt.legend()
pylab.tight_layout()
plt.savefig("evaluecoverage.png")
plt.savefig("evaluecoverage.pdf")
plt.show()
