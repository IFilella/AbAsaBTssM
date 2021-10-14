import fastaf
import numpy
import pandas as pd
import re
import os
import glob
import search
import subprocess

def do_proteinCheck(taxid,blastdb,protname,queries,ethr,lthr):
    #Look for protein homologous from a query list 'queries'
    for query in queries:
        out = query.split("/")[-1].split(".")[0]
        search.blastp(db=blastdb,query=query,out="%s/%s.%s.blasted"%(taxid,protname,out))
    catcommand = "cat %s/%s.*.blasted > %s/%s.blasted"%(taxid,protname,taxid,protname)
    os.system(catcommand)

    #Filter by evalue and lenght and convert to fasta
    search.filter_blastp_search("%s/%s.blasted"%(taxid,protname),"%s/%s.f.blasted"%(taxid,protname),evalue=ethr,lenght=lthr)
    search.get_fasta_from_blasted("%s/%s.f.blasted"%(taxid,protname),"%s/%s.f.fa"%(taxid,protname))
    output = subprocess.run("wc %s/%s.f.fa"%(taxid,protname),shell=True,capture_output=True)
    aux = str(output.stdout).split("\'")[1].split()[0]
    if int(aux) > 0: presence = "1"
    else: presence = "0"
    return presence

"""
###Add Organism and TaxID to each sequence title###

#Input data
taxiddb = "/work/ifilella/uniref90/uniref90.taxlist"
orgdb = "/work/ifilella/uniref90/names.dmp"
aliMctd = fastaf.Alignment("TssMctd/TssMctd.aln")

#Output data
f = open("TssMctd/TssMctd.complete.aln","w")

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
###Find GXXGXXXG pattern on non Acinetobacter TssMctd###

#Input data
alictd = fastaf.Alignment(aliname="TssMctd/TssMctd.complete.aln")
acinetotaxs = numpy.genfromtxt("data/acinetobacters.txids",dtype=str)
ettssj = 1e-20
lttssj = [100,500]
tssjs  = glob.glob("data/TssJ/*.TssJ.fa")
ettssb = 1e-20
lttssb = [100,250]
tssbs = glob.glob("data/TssB/*.TssB.fa")
queryAsaB = "data/AsaB/Ab.AsaB.fa"
genomes = "genomesTssMctd/" 

#Output data
f = open("TssMctd/GXXXGXXXG.noAcineto.txt","w")
f.write("GxxxGxxxG;key;organism;taxid;B;J;AsaB\n")

for i,name in enumerate(alictd.names):
    print("-----------------------------------------------------------------")
    print(i,name)
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
            continue
    #Filter by single genome association
    taxids = search.get_subTaxIDs(taxid)
    if len(taxids) > 1:
        # Not a leave TaxID can't recover a genome
        presenceJ = "Error5"
        presenceB = "Error5"
    else:
        if not os.path.isdir("%s%s"%(genomes,taxid)):
            os.system("mkdir %s%s"%(genomes,taxid))
        if os.path.isfile("%s%s/%s.fasta"%(genomes,taxid,taxid))==True:
            pass
        else:
            search.get_genome(taxid,"%s%s/%s.fasta"%(genomes,taxid,taxid))
        if os.stat("%s%s/%s.fasta"%(genomes,taxid,taxid)).st_size == 0:
            #Empty genome
            presenceJ = "Error6"
            presenceB = "Error6"
        else:
            search.make_blastdb("%s%s/%s.fasta"%(genomes,taxid,taxid),"%s%s/%s"%(genomes,taxid,taxid))
            blastdb = "%s%s/%s"%(genomes,taxid,taxid)
            #Look for TssJ homologs
            presenceJ = do_proteinCheck(taxid=genomes+taxid,blastdb,"TssJ",tssjs,ettssj,lttssj)
            presenceB = do_proteinCheck(taxid=genomes+taxid,blastdb,"TssB",tssbs,ettssb,lttssb)

    #Found the pattern on the sequence
    patterns1 = re.findall('G..G..G',unaliseq[-100:])
    patterns2 = re.findall('G...G...G',unaliseq[-100:])
    patterns3 = re.findall('G..G...G',unaliseq[-100:])
    patterns4 = re.findall('G...G..G',unaliseq[-100:])
    patterns = patterns1 + patterns2 + patterns3 + patterns4
    if len(patterns)>0:
        print("%s;%s;%s;%s;%s;%s\n"%('_'.join(patterns),key,organism,taxid,presenceB,presenceJ))
        f.write("%s;%s;%s;%s;%s;%s\n"%('_'.join(patterns),key,organism,taxid,presenceB,presenceJ))
    else:
        print("NO;%s;%s;%s;%s;%s\n"%(key,organism,taxid,presenceB,presenceJ))
        f.write("NO;%s;%s;%s;%s;%s\n"%(key,organism,taxid,presenceB,presenceJ))
f.close()
