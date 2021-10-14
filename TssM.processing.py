import glob
import fastaf
import search
import subprocess
import math
import os

#Input Data
taxids = glob.glob("genomesTssM/*")
orgdb = "/work/ifilella/uniref90/names.dmp"
genestinf = 500
genestsup = 8000
etTssM = 1e-20
coveragetTssM = 0.6
etAsaB = 1e-20
coveragetAsaB = 0.6
etTssJ = 1e-20
minltTssJ = 100
maxltTssJ = 500
etTssB = 1e-20
coveragetTssB = 0.6
etTssK = 1e-20
coveragetTssK = 0.6
queryAsaB = "data/AsaB/Ab.AsaB.fa"
queryTssM = "data/TssM/Ab.TssM.fa"
queryTssB = "data/TssB/Ab.TssB.fa"
queryTssK = "data/TssK/Ab.TssK.fa"
queriesTssJ = glob.glob("data/TssJ/*.TssJ.fa")

#Output Data
fanalysis = open("TssM/TssM.list.txt","w")
fanalysis.write("TaxID,Organism,TssM,AsaB,TssJ,TssB,TssK\n")
ftssm = open("TssM/Acinetobacter.TssM.fa","w")
fasab = open("TssM/Acinetobacter.AsaB.fa","w")
ftssb = open("TssM/Acinetobacter.TssB.fa","w")
ftssk = open("TssM/Acinetobacter.TssK.fa","w")


def get_protein(query,taxid,protname,eth,cth,fout,organism):
    tax = taxid.split("/")[1]
    queryseq = fastaf.fastaf(query).homolseqs[0].seq
    #Look for protein homologs#
    blasted = "%s/%s.Ab.blasted"%(taxid,protname)
    search.blastp(db=blastdb,query=query,out=blasted)
    #Apply a first filter based on e-value
    blastedf = "%s/%s.Ab.f.blasted"%(taxid,protname)
    search.filter_blastp_search(blasted,blastedf,evalue=eth)
    #Check that after the first filter there still are protein homologs and apply second coverage filter
    aux = search.get_lines(blastedf)
    if int(aux) > 0:
        search.filter_blastp_bycoverage(queryseq,blastedf,cth,blastedf)
        aux = search.get_lines(blastedf)
        if int(aux) > 0:
            fafile = "%s/%s.Ab.f.fa"%(taxid,protname)
            search.get_fasta_from_blasted(blastedf,fafile)
            presence = "1"
            #Get the homolog with the lowest evalue
            homols = fastaf.fastaf(fafile,searchtool="blastp").homolseqs
            eaux = math.inf
            prothomol = None
            for homol in homols:
                if homol.evalue < eaux:
                    eaux = homol.evalue
                    prothomol = homol.seq
            fout.write(">%s_%s\n"%(organism.replace(" ","_"),tax))
            fout.write("%s\n"%prothomol)
        else:
            presence = "0"
    else:
        presence = "0"
    return presence


count = 0
totalgenes = []
print(len(taxids))
for taxid in taxids:
    tax = taxid.split("/")[1]
    taxgenome = fastaf.fastaf("%s/%s.fasta"%(taxid,tax))
    taxgenes = len(taxgenome.homolseqs)
    totalgenes.append(taxgenes)
    if taxgenes > genestinf and taxgenes < genestsup:
        count+=1

        #Get the OrganismID
        output = subprocess.run("awk 'BEGIN {FS=\"\t\"} {if ($1==\"%s\"){print $3}}' %s"%(tax,orgdb), shell=True,capture_output=True)
        if output.returncode == 0:
            organism = str(output.stdout).split("\'")[1].split("\\")[0]
        else:
            #Organism not annotated
            organism="Error4"
        blastdb = "%s/%s"%(taxid,tax)

        presenceM = get_protein(query=queryTssM,taxid=taxid,protname="TssM",eth=etTssM,cth=coveragetTssM,fout=ftssm,organism=organism)
        presenceA = get_protein(query=queryAsaB,taxid=taxid,protname="AsaB",eth=etAsaB,cth=coveragetAsaB,fout=fasab,organism=organism)
        presenceB = get_protein(query=queryTssB,taxid=taxid,protname="TssB",eth=etTssB,cth=coveragetTssB,fout=ftssb,organism=organism)
        presenceK = get_protein(query=queryTssK,taxid=taxid,protname="TssK",eth=etTssK,cth=coveragetTssK,fout=ftssk,organism=organism)        

        ##Look for TssJ homologs##
        for TssJ in queriesTssJ:
            out = TssJ.split("/")[-1].split(".")[0]
            blasted = "%s/TssJ.%s.blasted"%(taxid,out)
            search.blastp(db=blastdb,query=TssJ,out=blasted)
        blasted = "%s/TssJ.blasted"%taxid
        catcommand = "cat %s/TssJ.*.blasted > %s"%(taxid,blasted)
        os.system(catcommand)
        #Filter by evalue and lenght and convert to fasta TssJ hits
        blastedf = "%s/TssJ.f.blasted"%taxid
        search.filter_blastp_search(blasted,blastedf,evalue=etTssJ,lenght=[minltTssJ,maxltTssJ])
        fafile = "%s/TssJ.f.fa"%taxid
        search.get_fasta_from_blasted(blastedf,fafile)
        aux = search.get_lines(blastedf)
        if int(aux) > 0 : presenceJ = 1
        else: presenceJ = 0
        
        print(count,tax,organism,presenceM,presenceA,presenceJ,presenceB,presenceK)
        fanalysis.write("%s;%s;%s;%s;%s;%s;%s\n"%(tax,organism,presenceM,presenceA,presenceJ,presenceB,presenceK))

print(count)
fanalysis.close()
ftssm.close()
fasab.close()
ftssb.close()
ftssk.close()
