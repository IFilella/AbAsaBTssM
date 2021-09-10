import glob
import numpy
import search
import fastaf
import subprocess
from time import sleep
import os

taxiddb = "/work/ifilella/uniref90/uniref90.taxlist"
orgdb = "/work/ifilella/uniref90/names.dmp"
tssms = glob.glob("data/TssM/*.TssM.fa")
tssjs = glob.glob("data/TssJ/*.TssJ.fa")
ettssm = 1e-40
tssmmin = 700
tssmmax = 1500
ettssj = 1e-25
tssjmin = 200
tssjmax = 700


f = open('AsaB/AsaB.txt',"w")
f.write("Homolog;Evalue;Coverage;TaxID;Organism;M;J\n")
titles = []
coverages = []
organisms = []
taxids = []
presenceJs = []
presenceMs = []
evalues = []

# Get Taxid,Organism,evalues and coverage for each homolog
homologs = "AsaB/Ab.AsaB.fa"
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
                else:
                    search.make_blastdb("genomesAsaB/%s/%s.fasta"%(taxid,taxid),"genomesAsaB/%s/%s"%(taxid,taxid))
                    blastdb = "genomesAsaB/%s/%s"%(taxid,taxid)
                
                    #Look for TssM homologous from a query list 'tssms'
                    for tssm in tssms:
                        out = tssm.split("/")[-1].split(".")[0]
                        search.blastp(db=blastdb,query=tssm,out="genomesAsaB/%s/TssM.%s.blasted"%(taxid,out))
                    catcommand = "cat genomesAsaB/%s/TssM.*.blasted > genomesAsaB/%s/TssM.blasted"%(taxid,taxid)
                    os.system(catcommand)
                    
                    #Filter by evalue and lenght and convert to fasta TssM hits
                    search.filter_blastp_search("genomesAsaB/%s/TssM.blasted"%taxid,"genomesAsaB/%s/TssM.f.blasted"%taxid,evalue=ettssm,lenght=[tssmmin,tssmmax])
                    search.get_fasta_from_blasted("genomesAsaB/%s/TssM.f.blasted"%taxid,"genomesAsaB/%s/TssM.f.fa"%taxid)
                    output = subprocess.run("wc genomesAsaB/%s/TssM.f.fa"%taxid,shell=True,capture_output=True)
                    aux = str(output.stdout).split("\'")[1].split()[0]
                    if int(aux) > 0: presenceM = "1"
                    else: presenceM = "0"

                    #Look for TssJ homologous from a query list 'tssjs'
                    for tssj in tssjs:
                        out = tssj.split("/")[-1].split(".")[0]
                        search.blastp(db=blastdb,query=tssj,out="genomesAsaB/%s/TssJ.%s.blasted"%(taxid,out))
                    catcommand = "cat genomesAsaB/%s/TssJ.*.blasted > genomesAsaB/%s/TssJ.blasted"%(taxid,taxid)
                    os.system(catcommand)
                    
                    #Filter by evalue and lenght and convert to fasta TssJ hits
                    search.filter_blastp_search("genomesAsaB/%s/TssJ.blasted"%taxid,"genomesAsaB/%s/TssJ.f.blasted"%taxid,evalue=ettssj,lenght=[tssjmin,tssjmax])
                    search.get_fasta_from_blasted("genomesAsaB/%s/TssJ.f.blasted"%taxid,"genomesAsaB/%s/TssJ.f.fa"%taxid)
                    output = subprocess.run("wc genomesAsaB/%s/TssJ.f.fa"%taxid,shell=True,capture_output=True)
                    aux = str(output.stdout).split("\'")[1].split()[0]
                    if int(aux) > 0: presenceJ = "1"
                    else: presenceJ = "0"
                    
        else:
            presenceM=taxid
            presenceJ=taxid
        presenceMs.append(presenceM)
        presenceJs.append(presenceJ)
        print("---------------------------------")
        print(tit,taxid,organism,presenceM,presenceJ)
        print("---------------------------------")

print(len(titles),len(evalues),len(coverages),len(taxids),len(organisms),len(presenceMs),len(presenceJs))
for i,title in enumerate(titles):
    f.write("%s;%s;%.3f;%s;%s;%s;%s\n"%(titles[i],evalues[i],coverages[i],taxids[i],organisms[i],presenceMs[i],presenceJs[i]))

