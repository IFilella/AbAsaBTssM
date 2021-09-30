import glob
import numpy
import search
import fastaf
import subprocess
from time import sleep
import os

def do_proteinCheck(protname,queries,ethr,lthr):
    #Look for protein homologous from a query list 'queries'
    for query in queries:
        out = query.split("/")[-1].split(".")[0]
        search.blastp(db=blastdb,query=query,out="genomesAsaB/%s/%s.%s.blasted"%(taxid,protname,out))
    catcommand = "cat genomesAsaB/%s/%s.*.blasted > genomesAsaB/%s/%s.blasted"%(taxid,protname,taxid,protname)
    os.system(catcommand)
    
    #Filter by evalue and lenght and convert to fasta
    search.filter_blastp_search("genomesAsaB/%s/%s.blasted"%(taxid,protname),"genomesAsaB/%s/%s.f.blasted"%(taxid,protname),evalue=ethr,lenght=lthr)
    search.get_fasta_from_blasted("genomesAsaB/%s/%s.f.blasted"%(taxid,protname),"genomesAsaB/%s/%s.f.fa"%(taxid,protname))
    output = subprocess.run("wc genomesAsaB/%s/%s.f.fa"%(taxid,protname),shell=True,capture_output=True)
    aux = str(output.stdout).split("\'")[1].split()[0]
    if int(aux) > 0: presence = "1"
    else: presence = "0"
    return presence

#Input Data
taxiddb = "/work/ifilella/uniref90/uniref90.taxlist"
orgdb = "/work/ifilella/uniref90/names.dmp"
tssms = glob.glob("data/TssM/*.TssM.fa")
tssjs = glob.glob("data/TssJ/*.TssJ.fa")
tssbs = glob.glob("data/TssB/*.TssB.fa")
tssks = glob.glob("data/TssK/*.TssK.fa")
ettssm = 1e-40
lttssm = [700,1500]
ettssj = 1e-30
lttssj = [200,700]
ettssb = 1e-30
lttssb = [100,250]
ettssk = 1e-30
lttssk = [300,650]
homologs = "AsaB/Ab.AsaB.fa"

#Output Data
f = open('AsaB/AsaB.2.txt',"w")
f.write("Homolog;Evalue;Coverage;TaxID;Organism;M;J;B;K\n")
titles = []
coverages = []
organisms = []
taxids = []
presenceMs = []
presenceJs = []
presenceBs = []
presenceKs = []
evalues = []

# Get Taxid,Organism,evalues and coverage for each homolog
homols = fastaf.fastaf(homologs)
lenquery = homols.get_query().length
query =  homols.get_query()
for i,homol in enumerate(homols.homolseqs):
    if "UniRef90" in homol.title:
        print(i,homol.key)
        tit=homol.title.split()[0]
        titles.append(tit)
        evalues.append(homol.title.split()[1])
        aux = fastaf.get_pairAlignment(seq1=query.seq,seq2=homol.seq,alimode='muscle')
        coverages.append(aux[1])
        output = subprocess.run("wget https://www.uniprot.org/uniprot/%s"%homol.key,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
        if output.returncode == 1:
            search.rm_getw(homol.key)
            #Impossible to recover the sequence from UNIPROT
            taxid="Error1"
            organism="Error1"
        else:
            output = subprocess.run("grep 'var isObsolete = true' %s"%homol.key, shell=True,capture_output=True)
            if output.returncode == 0:
                search.rm_getw(homol.key)
                #UNIPROT identifier obsolete
                taxid="Error2"
                organism="Error2"
            else:
                output = subprocess.run("grep 'OX=' %s"%(homol.key), shell=True,capture_output=True)
                if output.returncode == 0:
                    taxid = str(output.stdout).split("OX=")[1].split()[0]
                    organism = str(output.stdout).split("OS=")[1].split("OX=")[0]
                    search.rm_getw(homol.key)
                else:
                    output = subprocess.run("grep %s %s"%(homol.key,taxiddb), shell=True,capture_output=True)
                    if output.returncode == 0:
                        taxid = str(output.stdout).split("\\t")[1].split("\\")[0]
                        output = subprocess.run("awk 'BEGIN {FS=\"\t\"} {if ($1==\"%s\"){print $3}}' %s"%(taxid,orgdb), shell=True,capture_output=True)
                        if output.returncode == 0:
                            organism = str(output.stdout).split("\'")[1].split("\\")[0]
                            search.rm_getw(homol.key)
                        else:
                            #Organism not annotated
                            search.rm_getw(homol.key)
                            organism="Error4"
                    else:
                        #TaxID not annotated
                        search.rm_getw(homol.key)
                        taxid="Error3"
        organisms.append(organism)
        taxids.append(taxid)
        
        if "Error" not in taxid and "Error" not in organism:
            #Check if the taxid is a leave
            aux = search.get_subTaxIDs(taxid)
            if len(aux) > 1:
                # Not a leave TaxID can't recover a genome
                presenceM="Error5"
                presenceJ="Error5"
                presenceB="Error5"
                presenceK="Error5"
            else: 
                #Download the genome for the given taxid
                if not os.path.isdir("genomesAsaB/%s"%taxid):
                    os.system("mkdir genomesAsaB/%s"%taxid)
                if os.path.isfile("genomesAsaB/%s/%s.fasta"%(taxid,taxid))==True:
                    pass
                else:
                    search.get_genome(taxid,"genomesAsaB/%s/%s.fasta"%(taxid,taxid))
                if os.stat("genomesAsaB/%s/%s.fasta"%(taxid,taxid)).st_size == 0:
                    print('There are no sequences in the genome of the given taxid')
                    #Empty genome
                    presenceM="Error6"
                    presenceJ="Error6"
                    presenceB="Error5"
                    presenceK="Error5"
                else:
                    search.make_blastdb("genomesAsaB/%s/%s.fasta"%(taxid,taxid),"genomesAsaB/%s/%s"%(taxid,taxid))
                    blastdb = "genomesAsaB/%s/%s"%(taxid,taxid)
                    
                    presenceM = do_proteinCheck("TssM",tssms,ettssm,lttssm)
                    presenceJ = do_proteinCheck("TssJ",tssjs,ettssj,lttssj)         
                    presenceB = do_proteinCheck("TssB",tssbs,ettssb,lttssb)
                    presenceK = do_proteinCheck("TssK",tssks,ettssk,lttssk)

                    
        else:
            presenceM=taxid
            presenceJ=taxid
            presenceB=taxid
            presenceK=taxid
        presenceMs.append(presenceM)
        presenceJs.append(presenceJ)
        presenceBs.append(presenceB)
        presenceKs.append(presenceK)
        print("---------------------------------")
        print(tit,taxid,organism,presenceM,presenceJ,presenceB,presenceK)
        print("---------------------------------")

print(len(titles),len(evalues),len(coverages),len(taxids),len(organisms),len(presenceMs),len(presenceJs),len(presenceBs),len(presenceKs))
for i,title in enumerate(titles):
    f.write("%s;%s;%.3f;%s;%s;%s;%s;%s;%s\n"%(titles[i],evalues[i],coverages[i],taxids[i],organisms[i],presenceMs[i],presenceJs[i],presenceBs[i],presenceKs[i]))

