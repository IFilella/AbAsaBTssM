#!/usr/bin/env python3
# -*- coding: UTF8 -*-

import numpy as np
import os
import copy
from collections import OrderedDict
import hashlib
import random
from Bio.SubsMat import MatrixInfo as matlist
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import collections
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pylab
import subprocess
import progressbar
from time import sleep

class homolseq(object):
    """
    Init the object homologous sequence from a fasta file
    """
    def __init__(self,seq,title,searchtool,fulltitle,table):
        self.seq=seq
        self.title=title
        self.database=None
        self.key=None
        self.length=len(seq)
        self.score=None
        self.evalue=None
        self.query=None
        if "UniRef90" in title and fulltitle == False:
            self.key=self.title.split("UniRef90_")[1].split("/")[0].split()[0]
            self.database="uniref90"
            self.query=False
        elif "SGB" in title and fulltitle == False:
            self.key=self.title.split(":")[0]
            self.database="SGB"
            self.query=False
        elif "+" == title[-1] or "-"== title[-1] and fulltitle == False:
            self.key=title.split()[0]
            self.database="raptorX"
            self.query=False
        elif fulltitle == False:
            self.key=title.split()[0]
            self.database="query"
            self.query=True
        else:
            self.database="query"
            self.query=True
        if searchtool=="jackhmmer":
            if self.query:
                self.evalue=0
                self.score=0
            else:
                command="grep %s %s"%(self.key,table)
                out=os.popen(command).read()
                self.evalue=float(out.split()[4])
                self.score=float(out.split()[5])
        elif searchtool=="blastp":
            try:
                self.evalue=float(self.title.split()[1])
            except:
                self.evalue=0
            try:
                self.score=float(self.title.split()[2])
            except:
                self.evalue=0
        elif searchtool == None:
            self.evalue = None
            self.score = None
            

class Alignment(object):
    """
    Init the object alignment 
    """
    def __init__(self,aliname=None,dictionary=None,uniqueid=False):
        if dictionary is not None:
            self.ali = dictionary
        else:
            self.ali = OrderedDict()
        if aliname is not None:
            self.aliname=aliname
            self.read_alifile(self.aliname,uniqueid)
        else:
            self.aliname=None
        if self.ali:
            self.names = list(self.ali.keys())
            self.seqs = list(self.ali.values())
        else:
            self.names = []
            self.seqs = []
        if self.seqs != []:
            self.length = len(self.seqs[0])
            
    def get_trimmed_ali(self,selection,aliname=None):
        "selection must be a list with two integers (first ali position is 1)"
        trimmedAli=copy.copy(self.ali)
        for seq in trimmedAli:
            trimmedAli[seq]=trimmedAli[seq][selection[0]-1:selection[1]]
        return Alignment(aliname=aliname,dictionary=trimmedAli)
    
    def get_lenfilt_ali(self,):
        trimmedAli=copy.copy(self.ali)
        for seq in self.ali:
            print(seq)
        return 0

    def get_clean_EvalMSA_ali(self,Evalout,aliname=None,keep=None):
        cleanAli=copy.copy(self.ali)
        fclean=open(Evalout,"r")
        outliers=[]
        aux=False
        for l in fclean:
            if "ALIGNMENT ANALYSIS" in l:
                aux=False
            elif aux:
                line=l.split("(")[0]
                line=line.replace(" ","")
                line=line.replace("\n","")
                if line!="":
                    outliers.append(line)
            elif "Outliers:" in l:
                aux=True
                try: line=l.split()[1].split("(")[0]
                except: break
                line=line.replace(" ","")
                line=line.replace("\n","")
                if line!="":
                    outliers.append(line)
        outliers=set(outliers)
        for o in outliers:
            if "NODE" in o:
                o=o.split(":")[0].split("/")[0]
            elif "UniRef90" in o:
                o=o.split("_")[1].split("/")[0]
            if keep!=None and keep in o:
                continue
            else:
                del cleanAli[o]
        return Alignment(dictionary=cleanAli)

    def read_alifile(self,aliname,uniqueid=False):
        "uniqueid create uniqueid based on sequence"
        f=open(aliname,"r")
        sequence=""
        for i,l in enumerate(f):
            line=l.replace("\n","")
            try: line[0] == ">"
            except:
                print("IndexError: Line %s of %s is empty"%(i,aliname))
                continue
            if line[0] == ">":
                if len(sequence) > 0:
                    if not uniqueid:
                       self.ali[key]=sequence
                    else:
                       sequence_no_ali=sequence.replace("-","")
                       seqid = hashlib.md5(sequence_no_ali).hexdigest()
                       self.ali[key+"."+seqid]=sequence
                line=line.replace(">","")
                line=line.split()[0]
                if "/" in line:
                    line=os.path.dirname(line)
                if "UniRef90" in line:
                    line=line.split("UniRef90_")[1]
                elif "SGB" in line:
                    line=line.split(":")[0]
                key = line
                sequence=""
            else:
                sequence+=line
            if not uniqueid:
                self.ali[key]=sequence
            else:
                sequence_no_ali=sequence.replace("-","")
                seqid = hashlib.md5(sequence_no_ali).hexdigest()
                self.ali[key+"."+seqid]=sequence
        emptyseq = []
        for k in self.ali.keys():
            seq = self.ali[k]
            if "".join(set(seq)) == "-": emptyseq.append(k)
        if len(emptyseq) > 0: 
            print("%d empty sequences removed"%(len(emptyseq)))
        for k in emptyseq:
            del self.ali[k]
    
    def write_ali(self,output,upper=False):
        f=open(output,"w")
        for k in self.ali:
            f.write(">"+k+"\n")
            if upper:
                f.write(self.ali[k].upper()+"\n")
            else:
                f.write(self.ali[k]+"\n")
   
    def write_fasta(self,output,upper=False):
        f=open(output,"w")
        for k in self.ali:
            f.write(">"+k+"\n")
            if upper:
                f.write(self.ali[k].replace("-","").upper()+"\n")
            else:
                f.write(self.ali[k].replace("-","")+"\n")

    def do_patch_ids(self,Alignment2):
        """
        caveat: same size same sequence order
        """
        tmpali=OrderedDict()
        for n,k in enumerate(self.ali.keys()):
            newk=list(Alignment2.ali.keys())[n]
            tmpali[newk]=self.ali[k]
        self.ali=tmpali

    def rm_empty_columns(self):
        count=0
        todelete = []
        for i in range(self.length):
            delete = True
            for j,name in enumerate(self.names):
                if self.ali[name][i] != "-":
                    delete = False
                    break
            if delete:
                count+=1
                todelete.append(i)
        todelete.reverse()
        for i in todelete:
            for j,name in enumerate(self.names):
                self.ali[name] = self.ali[name][:i] + self.ali[name][i+1:]
        print("from %d sites %d were removed"%(self.length,count))
        self.length = len(self.seqs[0])
        self.seqs = self.ali.values()

    def get_concatenated_horizontally(self,Alignment2,spacer=""):
        """
        caveat: same size same ids order
        """
        alicopy=copy.copy(self.ali)
        for k in alicopy:
            seq2=Alignment2.ali[k]
            seq1=alicopy[k]
            alicopy[k]=seq1+spacer+seq2
        return Alignment(aliname=None,dictionary=alicopy)

    def replace_U_to_T(self,reverse=False):
        for seq in self.ali:
            if not reverse:
                self.ali[seq] = self.ali[seq].replace("U","T")
            if reverse:
                self.ali[seq] = self.ali[seq].replace("T","U")
 
    def __add__(self,Alignment2):
        alicopy=copy.copy(self.ali)
        print(Alignment2)
        alicopy.update(Alignment2.ali)
        return Alignment(aliname=None,dictionary=alicopy)


class fastafcont(object):
      """
      Init the object containing more than one fastaf that might come from different
      databases and diferent searchtools
      """
      def __init__(self,name):
          self.name=name
          self.searchtools=set()
          self.fastafs=[]
      
      def add_fastaf(self,fastaf):
          self.fastafs.append(fastaf)
          self.searchtools.add(fastaf.searchtool)
      
      def print_fasta(self,output):
          ffa=open(output,"w")
          if len(self.fastafs)!=0:
              for fa in self.fastafs:
                  for se in fa.homolseqs:
                      ffa.write(">%s\n"%se.title)
                      ffa.write("%s\n"%se.seq)

class fastaf(object): 
    """
    Init the object using a fasta file. It might be an homologous search done with
    jackhmmer or blast on uniref90 or SGB or it can be a regular fasta file.
    """
    def __init__(self,fastaname,searchtool=None,fulltitle=False,table=None):
        self.fastaname=fastaname
        self.searchtool=searchtool
        self.fulltitle=fulltitle
        self.homolseqs=[]
        self.table=table
        self.read_fasta()      
 
    def get_query(self):
        for seq in self.homolseqs:
            if seq.query:
                return seq
    
    def get_minseq_length(self):
        return int(self.get_query().length/2)

    def get_maxseq_length(self):
        return int(self.get_query().length*3/2)

    def read_fasta(self):
        homolseqs=[]
        f=open(self.fastaname,"r")
        count=0
        for l in f:
            l=l.replace("\n","")
            if l[0]==">":
                if count!=0:
                    homolseqs.append(homolseq(seq,title,self.searchtool,self.fulltitle,self.table))
                title=l.replace(">","")
                seq=""
                count+=1
            else:
                seq=seq+l
        homolseqs.append(homolseq(seq,title,self.searchtool,self.fulltitle,self.table))
        self.homolseqs=homolseqs
    
    def print_fasta(self,output,evalue=False):
        ffa=open(output,"w")
        for se in self.homolseqs:
            if evalue:
                ffa.write(">%s %.2e\n"%(se.title,se.evalue))
            else:
                ffa.write(">%s\n"%se.title)
            ffa.write("%s\n"%se.seq)
            

    def do_filtbyquerylength_fasta(self):
        filtindexes=[]
        for i,seq in enumerate(self.homolseqs):
            if self.get_minseq_length() < seq.length and seq.length < self.get_maxseq_length():
                pass
            else:
                filtindexes.append(i)
        filtindexes=set(filtindexes)
        filtindexes=list(filtindexes)
        filtindexes=sorted(filtindexes,reverse=True)
        for i in filtindexes:
            del self.homolseqs[i]

    def do_filted_fasta(self,ethr=None,lthr=[None,None]):
        filtindexes=[]
        for i,seq in enumerate(self.homolseqs):
            if ethr:
                if seq.evalue <= ethr or seq.evalue =="-":
                    pass
                else:
                    filtindexes.append(i)
            if lthr!=[None,None]:
                if lthr[0] < seq.length and seq.length < lthr[1]:
                    pass
                else:
                    filtindexes.append(i)
        filtindexes=set(filtindexes)
        filtindexes=list(filtindexes)
        filtindexes=sorted(filtindexes,reverse=True)
        for i in filtindexes:
            del self.homolseqs[i]

    def do_filtedbyquartile_fasta(self,quart=0.75):
        evalues=[]
        for i,seq in enumerate(self.homolseqs):
            evalues.append((i,seq.evalue))
        evalues=sorted(evalues, key=lambda tup: tup[1])
        totlen=len(evalues)
        evaluesfilt=evalues[int(quart*totlen):]
        filtindex=list(zip(*evaluesfilt))[0]
        filtindex=sorted(filtindex,reverse=True)
        for i in filtindex:
            del self.homolseqs[i]
        return 0
    
 
    def do_filtedbytax_fasta(self,tax):
        filtindexes=[]
        bar = progressbar.ProgressBar(maxval=len(self.homolseqs), \
            widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
        bar.start()
        for i,seq in enumerate(self.homolseqs):
            if seq.database=="uniref90":
                FNULL = open(os.devnull, 'w')
                output = subprocess.run("wget https://www.uniprot.org/uniprot/%s"%seq.key,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
                if output.returncode == 1:
                    subprocess.run("rm %s"%seq.key,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
                    bar.update(i+1)
                    sleep(0.001)
                    continue
                output = subprocess.run("grep 'var isObsolete = true' %s"%seq.key, shell=True,capture_output=True)
                if output.returncode == 0:
                    subprocess.run("rm %s"%seq.key,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
                    bar.update(i+1)
                    sleep(0.01)
                    continue
                output = subprocess.run("grep %s %s"%(tax,seq.key), shell=True,capture_output=True)
                if output.returncode == 1:
                    filtindexes.append(i)
                subprocess.run("rm %s"%seq.key,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
                bar.update(i+1)
                sleep(0.01)
        bar.finish()
        filtindexes=set(filtindexes)
        filtindexes=list(filtindexes)
        filtindexes=sorted(filtindexes,reverse=True)
        print("Number of no Proteobacterial seqs in homologs = " + str(len(filtindexes)))
        for i in filtindexes:
            del self.homolseqs[i]

    def get_homolseqs_keys(self):
        keys=[]
        for homolseq in self.homolseqs:
            keys.append(homolseq.key)
        return keys
    
    def get_homolseq(self,key):
        for homolseq in self.homolseqs:
            if homolseq.key == key:
                return homolseq.seq

    def get_title(self,key):
        for homolseq in self.homolseqs:
            if homolseq.key == key:
                return homolseq.title

    def do_shuffle_homolseqs(self):
        random.shuffle(self.homolseqs)

    def plot_IDheatmap(self,out,alimode="biopython",outali=None):
        dicc = collections.defaultdict(dict)
        if outali != None: fali = open(outali,"w")
        for i,seq1 in enumerate(self.homolseqs):
            for j,seq2 in enumerate(self.homolseqs):
                if i>=j: continue
                idy,ali = get_pairAlignment(seq1=seq1.seq,seq2=seq2.seq,gap_e=-1,gap_o=-5,alimode=alimode)
                dicc[seq1.key.split("|")[0]][seq2.key.split("|")[0]] = idy
                if outali != None:
                    if "Ab" in seq1.key:
                        print(ali)
                        fali.write(">"+seq1.key+"\n"+ali[0]+"\n"+">"+seq2.key+"\n"+ali[1]+"\n")
                        fali.write("-----------------------------------------------------\n")
        if outali != None: fali.close()
        df = pd.DataFrame(dicc)
        sns.heatmap(df,annot=True)
        plt.tight_layout()
        plt.savefig(out+'.png')
        pylab.show()

def get_pairAlignment(seq1,seq2,gap_e=-0.5,gap_o=-5,alimode="biopython"):
    matrix = matlist.blosum62
    if alimode == "biopython":
        alignment = pairwise2.align.globalds(seq1,seq2,matrix,gap_o,gap_e)
        alignment = np.asarray(alignment)
        scores = list(alignment[:,2].astype(float))
        max_score = max(scores)
        max_score_indx = scores.index(max(scores))
        alig1 = alignment[max_score_indx][0]
        alig2 = alignment[max_score_indx][1]
    elif alimode == "mafft" or alimode == "muscle" or "tcoffe":
        tmpfa = open("tmp.fa","w")
        tmpfa.write(">seq1"+"\n"+seq1+"\n"+">seq2"+"\n"+seq2)
        tmpfa.close()
        if alimode == "mafft":
            alicommand = "mafft tmp.fa > tmp.aln"
        elif alimode == "muscle":
            alicommand = "muscle -in tmp.fa -out tmp.aln"
        elif alimode == "tcoffe":
            alicommand = "t_coffee tmp.fa -output fasta -outfile tmp.aln"
        subprocess.run(alicommand,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
        os.system("rm tmp.fa")
        ali = Alignment("tmp.aln")
        alig1 = ali.ali["seq1"]
        alig2 = ali.ali["seq2"]
        os.system("rm tmp.aln")
    else:
        raise ValueError("Mode must be 'biopython', 'muscle' or 'mafft'")
    identity = 0
    for i in  range(len(alig1)):
        if (alig1[i] == alig2[i]) and alig1[i] != "-":
            identity +=1
    coverage = 0
    for i in range(len(alig1)):
        if alig1[i] !="-" and alig2[i] != "-":
            coverage +=1
    coverage  = float(coverage) / float(len(seq1))
    identity = float(identity)/len(alig1)
    return identity, coverage, [alig1,alig2]
