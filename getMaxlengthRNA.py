import sys
import re

import gzip
def readFile(infile):
    if infile.endswith((".gz","gzip")):
        fin = gzip.open(infile,'rt')
    else:
        fin = open(infile,'r')
    return fin
        
def writeFile(outfile):
    if outfile.endswith((".gz","gzip")):
        fout = gzip.open(outfile,'wt')
    else:
        fout = open(outfile,'w')
    return fout

from collections import OrderedDict
def parseString(mystr):
    d = OrderedDict()
    if mystr.endswith(";"):
        mystr = mystr[:-1]
    for i in mystr.strip().split(";"):
        tmp = i.strip().split()
        if tmp == "":
            continue
        d[tmp[0]] = tmp[1].replace("\"","")
    return d

genedict = {}
transdict = {}
fin = readFile(sys.argv[1])
for line in fin:
    if line.startswith("#"):
        continue
    tmp = line.strip().split("\t",8)
    if tmp[2] == "gene":
        d = parseString(tmp[8])
        if d["gene_id"] not in genedict:
            genedict[d["gene_id"]] = set()
        else:
            print(line)
            sys.exit("error1")
    if tmp[2] == "transcript":
        d = parseString(tmp[8])
        if d["gene_id"] not in genedict:
            print(line)
            sys.exit("error2")
        else:
            genedict[d["gene_id"]].add(d["transcript_id"])
        if d["transcript_id"] not in transdict:
            transdict[d["transcript_id"]] = {}
            transdict[d["transcript_id"]]["length"] = int(tmp[4]) - int(tmp[3]) + 1
        else:
            print(line)
            sys.exit("error3")
    if tmp[2] == "exon":
        d = parseString(tmp[8])
        if d["transcript_id"] not in transdict:
            transdict[d["transcript_id"]] = {}
            transdict[d["transcript_id"]]["exonlength"] = 0
        if "exonlength" not in transdict[d["transcript_id"]]:
            transdict[d["transcript_id"]]["exonlength"] = 0
        transdict[d["transcript_id"]]["exonlength"] += (int(tmp[4]) - int(tmp[3]) + 1)
fin.close()

outgene = set()
outtrans = set()
for k in genedict.keys():
    outgene.add(k)
    # print(len(genedict[k]))
    if len(genedict[k]) <= 3:
        for j in genedict[k]:
            outtrans.add(j)
    else:
        # top3 length
        translength = []
        for j in genedict[k]:
            translength.append([j,transdict[j]["length"]])
        translength.sort(key=lambda x:(x[1]),reverse=True)
        # top3 exon length
        transexonlength = []
        for j in genedict[k]:
            transexonlength.append([j,transdict[j]["length"]])
        transexonlength.sort(key=lambda x:(x[1]),reverse=True)
        outtrans.add(translength[0][0])
        outtrans.add(translength[1][0])
        outtrans.add(translength[2][0])
        outtrans.add(transexonlength[0][0])
        outtrans.add(transexonlength[1][0])
        outtrans.add(transexonlength[2][0])
# print(outtrans)

fin = readFile(sys.argv[1])
for line in fin:
    if line.startswith("#"):
        print(line,end="")
        continue
    tmp = line.strip().split("\t",8)
    if tmp[2] == "gene":
        d = parseString(tmp[8])
        if d["gene_id"] in outgene:
            print(line,end="")
    else:
        d = parseString(tmp[8])
        if d["transcript_id"] in outtrans:
            print(line,end="")
fin.close()