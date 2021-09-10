import os
import sys
import subprocess
from time import sleep
import fastaf

def get_genome(taxid,out):
    command = "esearch -db protein -query \"refseq[filter] AND txid%s[Organism]\" | efetch -format fasta > %s"%(taxid,out)
    os.system(command)

def make_blastdb(fasta,out):
    command = "makeblastdb -in %s -parse_seqids -title %s -dbtype prot -out %s"%(fasta,out,out)
    os.system(command)

def blastp(db,query,out,maxseqs=1000):
    command = "blastp -db %s -query %s -out %s -outfmt '6 sseqid sseq evalue score' -max_target_seqs %d"%(db,query,out,maxseqs)
    os.system(command)

def get_fasta_from_blasted(blasted,out):
    fin = open(blasted,"r")
    fout = open(out,"w")
    for line in fin:
        line = line.replace("\n","").split("\t")
        fout.write(">"+line[0].split("|")[1]+ " "+line[2]+" "+line[3]+"\n"+line[1]+"\n")
    fin.close()
    fout.close()

def filter_blastp_bycoverage(queryseq,blasted,coverageT,out):
    fin = open(blasted,"r")
    clines = []
    for line in fin:
        seq = line.split("\t")[1]
        aux = fastaf.get_pairAlignment(seq1=queryseq,seq2=seq,alimode='mafft')
        coverage = aux[1]
        if coverage > coverageT:
            clines.append(line.replace("\n",""))
    fin.close()
    fout = open(out,"w")
    for line in clines:
        fout.write(line+"\n")
    fout.close()

def filter_blastp_search(blasted,out,evalue=None,lenght=[None,None]):
    fin = open(blasted,"r")
    fout = open(out,"w")
    elines = []
    lenlines = []
    if evalue != None:
        for line in fin:
            if float(line.split("\t")[2])<=evalue:
                elines.append(line.replace("\n",""))
    else:
        for line in fin:
            elines.append(line.replace("\n",""))
    fin.close()
    fin = open(blasted,"r")
    if lenght[0] != None and lenght[1] != None:
        for line in fin:
            le = len(line.split("\t")[1].replace("-",""))
            if le >= int(lenght[0]) and le <= int(lenght[1]):
                lenlines.append(line.replace("\n",""))
    elif lenght[0] != None and lenght[1] == None:
        for line in fin:
            le = len(line.split("\t")[1].replace("-",""))
            if le >= int(lenght[0]):
                lenlines.append(line.replace("\n",""))
    elif lenght[0] == None and lenght[1] != None:
        for line in fin:
            le = len(line.split("\t")[1].replace("-",""))
            if le <= int(lenght[1]):
                lenlines.append(line.replace("\n",""))
    elif lenght[0] == None and lenght[1] == None:
        for line in fin:
            lenlines.append(line.replace("\n",""))
    fin.close()
    inter = list(set(elines).intersection(set(lenlines)))
    for line in inter:
        fout.write(line+"\n")
    fout.close()
    
def rm_getw(name):
    subprocess.run("rm %s*"%name,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
    sleep(0.001)

def get_lines(name):
    output = subprocess.run("wc %s"%name,shell=True,capture_output=True)
    aux = str(output.stdout).split("\'")[1].split()[0]
    sleep(0.001)
    return(aux)

def get_subTaxIDs(taxid,out=None):
    if out == None: filname = "tmp%s.txids"%taxid
    else: filname = "%s.txdis"%out
    subprocess.run("get_species_taxids.sh -t %s > %s"%(taxid,filname),shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL) 
    output = subprocess.run("wc %s"%filname,shell=True,capture_output=True)
    total = str(output.stdout).split("\'")[1].split()[0]
    f = open(filname,"r")
    taxids = []
    if int(total) >= 1:
        for line in f:
            taxids.append(line)
    if out == None: os.system('rm %s'%filname)
    return taxids
