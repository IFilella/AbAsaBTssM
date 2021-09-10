import search
import sys
import glob
import os

taxid = sys.argv[1]

##Download the genome for the given taxid##

#Check if the taxid have subtree TaxIDs
taxids = search.get_subTaxIDs(taxid)
if len(taxids) > 1:
    raise ValueError("TaxID with multiple subtree TaxIDs (Error5)")

#Download the genome if it wasn't done before
if not os.path.isdir("genomesTssM/%s"%taxid):
    os.system("mkdir genomesTssM/%s"%taxid)
if os.path.isfile("genomesTssM/%s/%s.fasta"%(taxid,taxid))==True:
    pass
else:
    search.get_genome(taxid,"genomesTssM/%s/%s.fasta"%(taxid,taxid))

#Check that the genome is not empty if not compute a Blast database
if os.stat("genomesTssM/%s/%s.fasta"%(taxid,taxid)).st_size == 0:
    os.system("rm -r genomesTssM/%s"%(taxid))
    raise ValueError ('There are no sequences in the genome of the given taxid (Error6)')
else:
    search.make_blastdb("genomesTssM/%s/%s.fasta"%(taxid,taxid),"genomesTssM/%s/%s"%(taxid,taxid))
    blastdb = "genomesTssM/%s/%s"%(taxid,taxid)

#Look for TssM homologs
#Look for AsaB homologs
