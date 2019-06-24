#!/usr/bin/env python3
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

opts, args = getopt.getopt(sys.argv[1:],"hf:m:c:l:",["folder=","mappingfile","readscount","length"])
input_file=""
def Usage():
    print ""
    print "This script can be used to estimate relative abundance of recovered MAGs."
    print "-h,  : Print help"
    print ""
    print "---------------------------------------------------------------------------------------------------------"
    print "Required options:"
    print "-f   : folder containing all MAGs"
    print "-m   : Mapping matrix eported from CLC or other mapping tools, which summarize the ID of contig/scaffold,"
    print "       mapped reads of given contig/scaffold, and length of this contig/scaffold. Example mapping matrix"
    print "       can be found in the folder of 'example-data'."
    print "-c   : which column (number) contains the mapped reads number."
    print "-l   : which column (number) contains the contigs/scaffolds length."
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
    elif op == "-l":
        columnL = value
        print "The ",columnL,"th column contains contig/scaffold length."
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



output1 = open("bins-average-coverage.txt","w") # file summarized all MAGs' average coverage
output2 = open("rel-abu-based-mapped-reads-of-bins.txt","w") #file summarized all MAGs' relative abundance (COVmag1/(COVmag1+...+COVmagn))
output1.write("Bin-id"+"\t"+"average-coverage"+"\t"+"Bin-size (bp)"+"\n")
output2.write("Bin-id"+"\t"+"Relative-abu"+"\n")


connections ={} #store the gene-ids and mapping results

#split mapping result into individual MAG
def getconnection(path,matrix):
    connections ={}
    tmp = open(matrix, "r")
    title = tmp.readline().strip()
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

getconnection(folder, mapping)
newfolder = "mapping-results-of-MAGs"
coveragefile = "bins-average-coverage.txt"

#get the average coverage of different MAGs
def Avecoverage(indivM,Counts,Length):
    col1 = int(Counts)-1
    col2 = int(Length)-1
    for root, dirs, files in os.walk(indivM):
        for file in files:
            tmp1 = np.loadtxt(os.path.join(root, file), skiprows=1, usecols=(col1, col2))
            num1 = np.sum(tmp1, axis=0)
            num1 = num1.tolist()
            rel = str(float(num1[0]) * 150 / float(num1[1]))
            output1.write(str(file).strip() + "\t" + str(rel).strip() + "\t" + str(num1[1]).strip() + "\n")
    output1.close()

#calculate the relative abundance: (average coverage of given MAG)/(total sum of all MAGs' coverage)
def Relabundance(BinCov):
    tmp2 = np.loadtxt(BinCov, skiprows=1, usecols=1)
    sum2 = np.sum(tmp2, axis=0)
    for line in open(BinCov,"r"):
        if "Bin-id" not in str(line):
            cont = str(line).strip().split("\t")
            relabu = round(float(cont[1])*100/float(sum2),3)
            output2.write(str(cont[0])+"\t"+str(relabu)+"\n")
        else:
            continue
    output2.close()


Avecoverage(newfolder,columnM,columnL)
Relabundance(coveragefile)


os.remove("bins-average-coverage.txt")

end = time.time()
print "Used time : ",str(end-start)
finish()
sys.exit()
