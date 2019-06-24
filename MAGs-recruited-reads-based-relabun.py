#!/usr/bin/env python2
###############################################################################
#    This script can be used to calculate the relative abundance of the       #
#    newly recovered MAGs.                                                    #
#    Example files for input can be found in the folder of "example-data"                                                                         #
#    This script was written and tested in python 2.7                         #
###############################################################################
__author__      = "Yulin Wang"
__maintainer__  = "Yulin Wang"
__email__       = "yulinwang605@gmail.com"
__status__      = "Development"

try:
    import numpy as np
except:
    print "Please install numpy."
try:
    from Bio import SeqIO
except:
    print "Please install biopython."

import os
import numpy as np
from Bio import SeqIO
import sys, getopt
import time

opts, args = getopt.getopt(sys.argv[1:],"hf:m:c:t:",["folder=","mappingfile","readscount","totalreads"])
input_file=""
def Usage():
    print ""
    print "This script can be used to estimate relative abundance of recovered MAGs (based on recruited reads number)."
    print "-h,  : Print help"
    print ""
    print "---------------------------------------------------------------------------------------------------------"
    print "Required options:"
    print "-f   : folder containing all MAGs"
    print "-m   : Mapping matrix exported from CLC or other mapping tools, which summarize the ID of contig/scaffold,"
    print "       mapped reads of given contig/scaffold, and length of this contig/scaffold. Example mapping matrix"
    print "       can be found in the folder of 'example-data'."
    print "-c   : which column (number) contains the mapped reads number."
    print "-t   : number of the total sequenced reads."
    print ""
def finish():
    print""
    print"Job finidhed!"
for op, value in opts:
    if op == "-f":
        folder = value
        print "The folder  is ", folder
    elif op == "-m":
        mapping = value
        print "The mapping matrix is ", mapping
    elif op == "-c":
        columnM = value
        print "The ",columnM,"th column contains number of mapped reads."
    elif op == "-t":
        Seqreads = value
        print "The number of total sequenced reads is: ",Seqreads
    elif op == "-h":
        Usage()
        sys.exit()

start = time.time()

if os.path.exists("mapping-results-of-MAGs"):
    for root, dirs, files in os.walk("mapping-results-of-MAGs"):
        for name in files:
            os.remove(os.path.join(root,name))
else:
    os.mkdir("mapping-results-of-MAGs")

output1 = open("bins-recruited-reads.txt","w") # file summarized all MAGs' average coverage
output2 = open("rel-abu-based-mapped-reads-count.txt","w") #file summarized all MAGs' relative abundance (COVmag1/(COVmag1+...+COVmagn))
output1.write("Bin-id"+"\t"+"average-coverage"+"\t"+"Bin-size (bp)"+"\n")
output2.write("Bin-id"+"\t"+"Relative-abun (vs. all sequenced reads)"+"\t"+"Relative-abun (vs. all reads mapped on assembled contigs"+"\n")

connections ={} #store the gene-ids and mapping results

#split mapping result into individual MAG
def getconnection(path,matrix,Counts):
    connections ={}
    tmp = open(matrix, "r")
    title = tmp.readline().strip()
    col1 = int(Counts) - 1
    totalMapped = np.loadtxt(matrix,skiprows=1, usecols=col1)
    totalMapped = np.sum(totalMapped,axis=0)
    for line in tmp:
        if "Name" not in str(line):
            cont = list(str(line).strip().split("\t"))
            connections[cont[0]]=str(cont[1]+'\t'+'\t'.join((cont[2:])))
        else:
            continue
    for root,dirs,files in os.walk(path):
        for file in files:
            output = open("mapping-results-of-MAGs/mapping-"+str(file)+".txt","w")
            output.write(str(title)+"\n")
            MAGfile = open(os.path.join(root,file),"r")
            for record in SeqIO.parse(MAGfile,"fasta"):
                id = str(record.id).strip()
                for key,value in connections.items():
                    if str(key)==str(id):
                        output.write(str(key)+"\t"+str(value)+"\n")
                    else:
                        continue
            output.close()
    return totalMapped

Mreads = getconnection(folder,mapping,columnM) #reads number that recruited by all the assembled contigs
print "Total number of reads mapped on the assembled contigs: ",Mreads

newfolder = "mapping-results-of-MAGs"
readscount = "bins-recruited-reads.txt"

#get the recruited reads of different MAGs
def RecruitR(indivM,Counts):
    col1 = int(Counts)-1
    for root, dirs, files in os.walk(indivM):
        for file in files:
            tmp1 = np.loadtxt(os.path.join(root, file), skiprows=1, usecols=col1)
            num1 = np.sum(tmp1, axis=0)
            output1.write(str(file).strip() + "\t" + str(float(num1)) + "\n")
    output1.close()

#calculate the relative abundance: (recrutied reads by given MAG/(total sequenced reads) or (recruited reads by given MAG)/(reads that mapped on all assembled contigs)
def Relabundance(BinCov):
    num = 0
    for line in open(BinCov,"r"):
        if "Bin-id" not in str(line):
            cont = str(line).strip().split("\t")
            num+=float(cont[1])
            relabu1 = round(float(cont[1])*100/float(Seqreads),3)
            relabu2 = round(float(cont[1])*100/float(Mreads),3)
            output2.write(str(cont[0])+"\t"+str(relabu1)+"\t"+str(relabu2)+"\n")
        else:
            continue
    unbinned1 = 100-round(num*100/float(str(Seqreads)),3)
    unbinned2 = 100-round(num*100/float(str(Mreads)),3)
    output2.write("unbinned"+"\t"+str(unbinned1)+"\t"+str(unbinned2)+"\n")
    output2.close()

RecruitR(newfolder,columnM)
Relabundance(readscount)

os.remove("bins-recruited-reads.txt")

end = time.time()
print "Used time : ",str(end-start)
finish()
sys.exit()