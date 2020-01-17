#!/usr/bin/python
# coding: UTF-8
import os


class Read:
    name = ''
    raw_string = ''
    pos_start_exon = []
    len_exon = []
    num_exon = 0
    start = 0
    end = 0
    score = 0
    strand = ''
    chrom = 0


def read(PATH):
    f = open(PATH)
    lines = f.readlines()
    f.close()
    reads = []
    for line in lines:
        ll = line.split("\t")
        R = Read()
        R.chrom = ll[0]
        if ll[0] == "chrX":
            R.chrom = "chr23"
        elif ll[0] == "chrY":
            R.chrom = "chr24"
        R.start = int(ll[1])
        R.end = int(ll[2])
        R.name = ll[3]
        R.score = round(float(ll[4]))
        R.strand = ll[5]
        R.cds_start = int(ll[6])
        R.cds_end = int(ll[7])
        R.itemRgb = ll[8]
        if(len(ll) > 9):
            R.num_exon = int(ll[9])
            R.size_exon = ll[10].rstrip(",").rstrip("").split(",")
            R.size_exon = list(map(int, R.size_exon))
            ll[11] = ll[11].rstrip("\n")
            R.pos_start_exon = ll[11].rstrip(",").split(",")
            R.pos_end_exon = []
            R.pos_start_exon = list(map(int, R.pos_start_exon))
            R.pos_start_exon = list(map(lambda x: x+R.start, R.pos_start_exon))
        for j in range(len(R.pos_start_exon)):
            try:
                R.pos_end_exon.append(
                    int(R.pos_start_exon[j])+int(R.size_exon[j]))
            except ValueError:
                R.pos_end_exon.append(0)
                print("illegal pos_end_exon")
        reads.append(R)
    return reads


def split_per_chrom(reads):
    reads_per_chrom = {}
    for read in reads:
        if read.chrom not in reads_per_chrom.keys():
            reads_per_chrom[read.chrom] = []
        reads_per_chrom[read.chrom].append(read)
    return reads_per_chrom


def execute_per_chrom(reads_per_chrom, func, **kwargs):
    # we don't execute at chrM, ambiguous chrom now
    output = {}
    for i in range(24):
        if("chr"+str(i+1) in reads_per_chrom):
            output["chr"+str(i+1)]=func(reads_per_chrom["chr"+str(i+1)], **kwargs)
    return output


def compare_per_chrom(reads1_per_chrom, reads2_per_chrom, func, **kwargs):
    # we don't execute at chrM, ambiguous chrom now
    output = {}
    for i in range(24):
        if("chr"+str(i+1) in reads1_per_chrom and
           "chr"+str(i+1) in reads2_per_chrom):
            output["chr"+str(i+1)]=func(reads1_per_chrom["chr"+str(i+1)],
                 reads2_per_chrom["chr"+str(i+1)], **kwargs)
    return output

def write(PATH, reads):
    f = open(PATH, "w")
    for read in reads:
        if read.chrom == "chr23":
            read.chrom = "chrX"
        elif read.chrom == "chr24":
            read.chrom = "chrY"
        read.pos_start_exon = ",".join(str(x) for x in read.pos_start_exon)
        read.size_exon = ",".join(str(x) for x in read.size_exon)
        f.write(read.chrom+"\t"+str(read.start)+"\t"+str(read.end)+"\t"+str(read.name)+"\t"
                + str(read.score)+"\t"+str(read.strand) +
                "\t"+read.cds_start+"\t"
                + read.cds_end+"\t"+"0"+"\t" +
                str(read.num_exon)+"\t"+read.size_exon+"\t"
                + read.pos_start_exon+"\n")
    f.close()
