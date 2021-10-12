import fastaf
import numpy
import pandas as pd
import re
import os
import glob
import search
import subprocess

def do_proteinCheck(protname,queries,ethr,lthr):
    #Look for protein homologous from a query list 'queries'
    for query in queries:
        out = query.split("/")[-1].split(".")[0]
        search.blastp(db=blastdb,query=query,out="genomesTssMctd/%s/%s.%s.blasted"%(taxid,protname,out))
    catcommand = "cat genomesTssMctd/%s/%s.*.blasted > genomesTssMctd/%s/%s.blasted"%(taxid,protname,taxid,protname)
    os.system(catcommand)

    #Filter by evalue and lenght and convert to fasta
    search.filter_blastp_search("genomesTssMctd/%s/%s.blasted"%(taxid,protname),"genomesTssMctd/%s/%s.f.blasted"%(taxid,protname),evalue=ethr,lenght=lthr)
    search.get_fasta_from_blasted("genomesTssMctd/%s/%s.f.blasted"%(taxid,protname),"genomesTssMctd/%s/%s.f.fa"%(taxid,protname))
    output = subprocess.run("wc genomesTssMctd/%s/%s.f.fa"%(taxid,protname),shell=True,capture_output=True)
    aux = str(output.stdout).split("\'")[1].split()[0]
    if int(aux) > 0: presence = "1"
    else: presence = "0"
    return presence

"""
###Add Organism and TaxID to each sequence title###

#Input data
taxiddb = "/work/ifilella/uniref90/uniref90.taxlist"
orgdb = "/work/ifilella/uniref90/names.dmp"
aliMctd = fastaf.Alignment("TssM/TssMctd.aln")

#Output data
f = open("TssM/TssMctd.complete.aln","w")

for i,name in enumerate(aliMctd.names):
    print("----------------")        
    if "TssM" in name:
        print(i+1,name)
        f.write(">%s\n"%name)
        f.write("%s\n"%aliMctd.seqs[i])
        continue
    output = subprocess.run("wget https://www.uniprot.org/uniprot/%s"%name,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
    if output.returncode == 1:
        search.rm_getw(name)
        #Impossible to recover the sequence from UNIPROT
        taxid="Error1"
        organism="Error1"
    else:
        output = subprocess.run("grep 'var isObsolete = true' %s"%name, shell=True,capture_output=True)
        if output.returncode == 0:
            search.rm_getw(name)
            #UNIPROT identifier obsolete
            taxid="Error2"
            organism="Error2"
        else:
            output = subprocess.run("grep 'OX=' %s"%(name), shell=True,capture_output=True)
            if output.returncode == 0:
                taxid = str(output.stdout).split("OX=")[1].split()[0]
                organism = str(output.stdout).split("OS=")[1].split("OX=")[0]
                search.rm_getw(name)
            else:
                output = subprocess.run("grep %s %s"%(name,taxiddb), shell=True,capture_output=True)
                if output.returncode == 0:
                    taxid = str(output.stdout).split("\\t")[1].split("\\")[0]
                    output = subprocess.run("awk 'BEGIN {FS=\"\t\"} {if ($1==\"%s\"){print $3}}' %s"%(taxid,orgdb), shell=True,capture_output=True)
                    if output.returncode == 0:
                        organism = str(output.stdout).split("\'")[1].split("\\")[0]
                        search.rm_getw(name)
                    else:
                        #Organism not annotated
                        search.rm_getw(name)
                        organism="Error4"
                else:
                    #TaxID not annotated
                    search.rm_getw(name)
                    taxid="Error3"
    print(i,name,organism,taxid)
    if "Error" not in organism and "Error" not in taxid:
        f.write(">%s|%s|%s\n"%(name,organism.replace(" ","_"),taxid))
        f.write("%s\n"%aliMctd.seqs[i])
f.close()
"""
"""
###Find GXXGXXXG pattern on non Acinetobacter TssMctd###

#Input data
alictd = fastaf.Alignment(aliname="TssM/TssMctd.complete.aln")
acinetotaxs = numpy.genfromtxt("TssM/acinetobacters.txids",dtype=str)
ettssj = 1e-30
lttssj = [200,700]
tssjs  = glob.glob("data/TssJ/*.TssJ.fa")

#Output data
total = 0
total1 = 0
total2 = 0
total3 = 0
totalplus = 0
totalAcineto = 0
f = open("TssM/GXXXGXXXG.txt","w")

for i,name in enumerate(alictd.names):
    unaliseq = alictd.seqs[i].replace("-","")
    aliseq = alictd.seqs[i]
    #Get uniprot identifier, organism and taxid
    if "TssM" not in name:
        try:
            taxid = name.split("|")[2]
        except: 
            taxid = ""
        organism = name.split("|")[1]
        key = name.split("|")[0]
    else:
        key = name.split("|")[1]
        organism = name.split("|")[0]
        taxid = ""
    #Get rid of the Acinetobacter sequences
    if taxid in acinetotaxs or 'Acineto' in name or 'Ab-TssM' in name:
            totalAcineto+=1
            continue 
    #Found the pattern on the sequence
    patterns1 = re.findall('G..G..G',unaliseq[-100:])
    patterns2 = re.findall('G...G...G',unaliseq[-100:])
    patterns3 = re.findall('G..G...G',unaliseq[-100:])
    patterns4 = re.findall('G...G..G',unaliseq[-100:])
    patterns = patterns1 + patterns2 + patterns3 + patterns4
    if len(patterns)>0:
        total+=1
        if len(patterns)==1: total1+=1
        if len(patterns)==2: total2+=1
        if len(patterns)==3: total3+=1
        if len(patterns)>3: totalplus+=1
        print(i,name,patterns)    
        
        #Download the genome associated to this specific taxid (if possible) to look for TssJ
        #Check if the taxid have subtree TaxIDs
        taxids = search.get_subTaxIDs(taxid)
        if len(taxids) > 1:
            # Not a leave TaxID can't recover a genome
            presenceJ = "Error5"
        else:
            if not os.path.isdir("genomesTssMctd/%s"%taxid):
                os.system("mkdir genomesTssMctd/%s"%taxid)
            if os.path.isfile("genomesTssMctd/%s/%s.fasta"%(taxid,taxid))==True:
                pass
            else:
                search.get_genome(taxid,"genomesTssMctd/%s/%s.fasta"%(taxid,taxid))
            if os.stat("genomesTssMctd/%s/%s.fasta"%(taxid,taxid)).st_size == 0:
                #Empty genome
                presenceJ = "Error6"
            else:
                search.make_blastdb("genomesTssMctd/%s/%s.fasta"%(taxid,taxid),"genomesTssMctd/%s/%s"%(taxid,taxid))
                blastdb = "genomesTssMctd/%s/%s"%(taxid,taxid)

                #Look for TssJ homologs
                presenceJ = do_proteinCheck("TssJ",tssjs,ettssj,lttssj)

        print("%s;%s;%s;%s;%s\n"%('_'.join(patterns),key,organism,taxid,presenceJ))
        f.write("%s;%s;%s;%s;%s\n"%('_'.join(patterns),key,organism,taxid,presenceJ))
f.close()
print(total,totalAcineto)
print(len(alictd.names),total,total1,total2,total3,totalplus)
"""
"""
###Analysis of the GxxxGxxxG output###

#Input data
gxfile = numpy.genfromtxt("TssM/GXXXGXXXG.txt",dtype=str,skip_header=1,delimiter=";")
print(len(gxfile))

#Output data
count = 0
tssjs = 0

for line in gxfile:
    print(line)
    if "Error" in line[-1]: continue
    else:
        count+=1
        if line[-1] == "0": tssjs+=1
print(count,tssjs)
"""
"""
###Find GXXGXXXG pattern on TssMctd of Acinetobacter species###

#Input data
Blist = numpy.genfromtxt("TssM/TssM.list.2.txt",skip_header=1,delimiter=";",dtype=str)
#alictd = fastaf.Alignment(aliname="TssM/Acinetobacter.TssM.2.fa")
#alictd = fastaf.Alignment(aliname="TssM/Baumannii.TssM.2.fa")
alictd = fastaf.Alignment(aliname="TssM/NOBaumannii.TssM.2.fa")

#Output data
count = 0
total = 0
total1 = 0
total2 = 0
total3 = 0
totalplus = 0
#f = open("TssM/GXXXGXXXG.Acineto.2.txt","w")
#f = open("TssM/GXXXGXXXG.Baumannii.2.txt","w")
f = open("TssM/GXXXGXXXG.NOBaumannii.2.txt","w")

for i,name in enumerate(alictd.names):
    unaliseq = alictd.seqs[i]
    taxid = name.split("_")[-1]
    organism = "_".join(name.split("_")[0:-2])
    indx = numpy.where(Blist==taxid)[0][0]
    presenceB = list(Blist[indx])[5]
    if presenceB != '1': continue
    count+=1
    #Found the pattern on the sequence
    patterns1 = re.findall('G..G..G',unaliseq[-100:])
    patterns2 = re.findall('G...G...G',unaliseq[-100:])
    patterns3 = re.findall('G..G...G',unaliseq[-100:])
    patterns4 = re.findall('G...G..G',unaliseq[-100:])
    patterns = patterns1 + patterns2 + patterns3 + patterns4
    if len(patterns)>0:
        total+=1
        print(i,name,patterns)
        if len(patterns)==1: total1+=1
        if len(patterns)==2: total2+=1
        if len(patterns)==3: total3+=1
        if len(patterns)>3: totalplus+=1
        print("%s;%s;%s\n"%('_'.join(patterns),organism,taxid))
        f.write("%s;%s;%s\n"%('_'.join(patterns),organism,taxid))
    else:
        f.write("NOTFOUND;%s;%s\n"%(organism,taxid))
print(len(alictd.names),count,total,total1,total2,total3,totalplus)
"""
"""
###Find GXXGXXXG pattern on TssMctd of Acinetobacter species###

#Input data
Blist = numpy.genfromtxt("TssM/TssM.list.2.txt",skip_header=1,delimiter=";",dtype=str)
#alictd = fastaf.Alignment(aliname="TssM/Acinetobacter.TssM.2.fa")
#alictd = fastaf.Alignment(aliname="TssM/Baumannii.TssM.2.fa")
alictd = fastaf.Alignment(aliname="TssM/NOBaumannii.TssM.2.fa")

#Output data
count = 0
total = 0
total1 = 0
total2 = 0
total3 = 0
totalplus = 0
#f = open("TssM/GXXXGXXXG.Acineto.2.txt","w")
#f = open("TssM/GXXXGXXXG.Baumannii.2.txt","w")
f = open("TssM/GXXXGXXXG.NOBaumannii.2.txt","w")

for i,name in enumerate(alictd.names):
    unaliseq = alictd.seqs[i]
    taxid = name.split("_")[-1]
    organism = "_".join(name.split("_")[0:-2])
    indx = numpy.where(Blist==taxid)[0][0]
    presenceB = list(Blist[indx])[5]
    if presenceB != '1': continue
    count+=1
    #Found the pattern on the sequence
    patterns1 = re.findall('G..G..G',unaliseq[-100:])
    patterns2 = re.findall('G...G...G',unaliseq[-100:])
    patterns3 = re.findall('G..G...G',unaliseq[-100:])
    patterns4 = re.findall('G...G..G',unaliseq[-100:])
    patterns = patterns1 + patterns2 + patterns3 + patterns4
    if len(patterns)>0:
        total+=1
        print(i,name,patterns)
        if len(patterns)==1: total1+=1
        if len(patterns)==2: total2+=1
        if len(patterns)==3: total3+=1
        if len(patterns)>3: totalplus+=1
        print("%s;%s;%s\n"%('_'.join(patterns),organism,taxid))
        f.write("%s;%s;%s\n"%('_'.join(patterns),organism,taxid))
    else:
        f.write("NOTFOUND;%s;%s\n"%(organism,taxid))
print(len(alictd.names),count,total,total1,total2,total3,totalplus)
"""
