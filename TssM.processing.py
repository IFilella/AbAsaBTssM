import glob
import fastaf
import search
import subprocess
import math
import os

taxids = glob.glob("genomesTssM/*")
genestinf = 500
genestsup = 8000
etTssM = 1e-25
coveragetTssM = 0.6
etAsaB = 1e-25
coveragetAsaB = 0.6
etTssJ = 1e-25
minltTssJ = 200
maxltTssJ = 700
queryAsaB = "/work/ifilella/AbAsaBTssM/AsaB/Ab.AsaB.fa"
queryTssM = "/work/ifilella/AbAsaBTssM/TssM/Ab.TssM.fa"
queriesTssJ = glob.glob("data/TssJ/*.TssJ.fa")
queryseqTssM = fastaf.fastaf(queryTssM).homolseqs[0].seq
queryseqAsaB = fastaf.fastaf(queryAsaB).homolseqs[0].seq
fanalysis = open("TssM/TssM.txt","w")
ftssm = open("TssM/Acinetobacters.TssM.fa","w")
fasab = open("TssM/Acinetobacters.AsaB.fa","w")
fanalysis.write("TaxID,Organism,TssM,AsaB,TssJ\n")
orgdb = "/work/ifilella/uniref90/names.dmp"

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

        ##Look for TssM homologs##
        blasted = "%s/TssM.Ab.blasted"%taxid
        search.blastp(db=blastdb,query=queryTssM,out=blasted)
        blastedf = "%s/TssM.Ab.f.blasted"%taxid
        search.filter_blastp_search(blasted,blastedf,evalue=etTssM)
        #Check that after the first filter there still are TssM homologs and apply second filter
        aux = search.get_lines(blastedf)
        if int(aux) > 0:
            search.filter_blastp_bycoverage(queryseqTssM,blastedf,coveragetTssM,blastedf)
            aux = search.get_lines(blastedf)
            if int(aux) > 0:
                fafile = "%s/TssM.Ab.f.fa"%taxid
                search.get_fasta_from_blasted(blastedf,fafile)
                presenceM = "1"
                #Get TssM homolog with lowest evalue
                homols = fastaf.fastaf(fafile,searchtool="blastp").homolseqs
                eaux = math.inf
                tssmhomol = None
                for homol in homols:
                    if homol.evalue < eaux:
                        eaux = homol.evalue
                        tssmhomol = homol.seq
                ftssm.write(">%s_%s\n"%(organism.replace(" ","_"),tax))
                ftssm.write("%s\n"%tssmhomol)
            else:
                presenceM = "0"
        else:
            presenceM = "0"

        ##Look for AsaB homologs##
        blasted = "%s/AsaB.Ab.blasted"%taxid
        search.blastp(db=blastdb,query=queryAsaB,out=blasted)
        blastedf = "%s/AsaB.Ab.f.blasted"%taxid
        search.filter_blastp_search(blasted,blastedf,evalue=etAsaB)
        #Check that after the first filter there still are TssM homologs and apply second filter
        aux = search.get_lines(blastedf)
        if int(aux) > 0:
            search.filter_blastp_bycoverage(queryseqAsaB,blastedf,coveragetAsaB,blastedf)
            aux = search.get_lines(blastedf)
            if int(aux) > 0:
                fafile = "%s/AsaB.Ab.f.fa"%taxid
                search.get_fasta_from_blasted(blastedf,fafile)
                presenceA = "1"
                #Get AsaB homolog with lowest evalue
                homols = fastaf.fastaf(fafile,searchtool="blastp").homolseqs
                eaux = math.inf
                asabhomol = None
                for homol in homols:
                    if homol.evalue < eaux:
                        eaux = homol.evalue
                        asabhomol = homol.seq
                fasab.write(">%s_%s\n"%(organism.replace(" ","_"),tax))
                fasab.write("%s\n"%asabhomol)
            else:
                presenceA = "0"
        else:
            presenceA = "0"

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
    
        print(count,tax,organism,presenceM,presenceA,presenceJ)
        fanalysis.write("%s;%s;%s;%s;%s\n"%(tax,organism,presenceM,presenceA,presenceJ))

print(count)
fanalysis.close()
ftssm.close()
fasab.close()
