#!/usr/bin/python
# coding: UTF-8

import os
import sys
import numpy as np
import bed_utils as bed
import matplotlib.pyplot as plt
import bed_utils as bed
from scipy.signal import argrelextrema
import pandas as pd
import seaborn as sns
import bisect
import argparse
from scipy import stats

def find_nearest_CAGEpeak(reads, cagereads):
    pTSS=[]
    mTSS=[]
    newreads = []
    f = open(args.input+".cagelog","a")
    for cageread in cagereads:
        if cageread.strand == "+":
            pTSS.append(cageread.cds_start)
        else:
            mTSS.append(cageread.cds_end)
    for read in reads:
        if read.strand == "+":
            minv = nearest_value(pTSS,read.start,args.tss_t)
            tss_mins.append(abs(abs(minv)-read.start))
            read.start = minv
        else:
            minv = nearest_value(mTSS,read.end,args.tss_t)
            tss_mins.append(abs(abs(minv)-read.end))
            read.end = minv
        if(minv>=0):
            newreads.append(read)
    return newreads

def find_nearest_splice_junction(reads,ss):
    newreads = []
    ss.sort()
    f = open(args.input+".sslog","a")
    for read_index,read in enumerate(reads):
        flag = 0
        for i in range(len(read.pos_start_exon)-1):
            minv = nearest_value(ss,read.pos_start_exon[i+1],args.ss_t)
            ss_mins.append(abs(abs(minv)-read.pos_start_exon[i+1]))
            if(minv>=0):
                read.pos_start_exon[i+1]=minv
            else:
                flag=1
        for i in range(len(read.pos_end_exon)-1):
            minv = nearest_value(ss,read.pos_end_exon[i],args.ss_t)
            ss_mins.append(abs(abs(minv)-read.pos_end_exon[i]))
            if(minv>=0):
                read.pos_end_exon[i]=minv
            else:
                flag=1
        if(read.strand == "+"):
            minv = nearest_value(ss,read.end,args.ss3end_t)
            ss3end_mins.append(abs(abs(minv)-read.end))
            if(minv>=0):
                read.end=minv
            else:
                flag=1
        else:
            minv = nearest_value(ss,read.start,args.ss3end_t)
            ss3end_mins.append(abs(abs(minv)-read.start))
            if(minv>=0):
                read.start=minv
            else:
                flag=1
        if(flag==0):
            newreads.append(read)
    return newreads

def nearest_value(ss,x,T):
    idx = bisect.bisect_left(ss,x)
    if(idx==0):
        if abs(ss[idx]-x) <= T:
            return ss[idx]
        else:
            return -1
    elif(idx==len(ss)):
        if abs(ss[idx-1]-x) <= T:
            return ss[idx-1]
        else:
            return -1
    if abs(ss[idx]-x)<abs(ss[idx-1]-x):
        minv = ss[idx]
    else:
        minv = ss[idx-1]
    if abs(minv-x) <= T:
        return minv
    else:
        return -1*minv

def create_sslist_from_genepred_file():
    f = open(args.ref)
    lines = f.readlines()
    ssperchrm = {}
    for line in lines:
        line = line.split("\t")
        chrm = line[2]
        if chrm[3:] == "X":
            chrm = "chr23"
        elif chrm[3:] == "Y":
            chrm = "chr24"
        elif not chrm[3:].isdigit():
            continue
        ss = list(map(int,line[9].split(",")[:-1])) + list(map(int,line[10].split(",")[:-1]))
        if chrm not in ssperchrm:
            ssperchrm[chrm] = ss
        else:
            ssperchrm[chrm] += ss
    return ssperchrm

def create_sslist_from_gtf_file():
    # currently accepts only chr, start, end format
    f = open(args.ref+".preprocessed")
    lines = f.readlines()
    ssperchrm = {}
    ssperchrm = {}
    for line in lines:
        line = line.split(" ")
        if(line[3]!="exon\n"):
            continue
        chrm = line[0]
        if chrm[3:] == "X":
            chrm = "chr23"
        elif chrm[3:] == "Y":
            chrm = "chr24"
        elif not chrm[3:].isdigit():
            continue
        ss = [int(line[1]),int(line[2])]
        if chrm not in ssperchrm:
            ssperchrm[chrm] = ss
        else:
            ssperchrm[chrm] += ss
    return ssperchrm


def input_to_ssdict():
    # identify format
    f = open(args.ref)
    line = f.readline()
    while line[0]=="#":
        line = f.readline()
    row = line.split(" ")
    if row[2][:3] == "chr":
        # genepred format
        ssperchrm = create_sslist_from_genepred_file()
    elif row[0][:3] == "chr":
        # gtf format
        if os.path.exists(args.ref+".preprocessed") == False:
            os.system("cat "+args.ref+" | awk '{print $1,$4,$5,$3}' > "+args.ref+".preprocessed")
        ssperchrm = create_sslist_from_gtf_file()
    else:
        raise "illegal_input"
    return ssperchrm

def write_corrected_ss(reads, **kwargs):
    OUTPUT_PATH = kwargs["OUTPUT_PATH"]
    f = open(OUTPUT_PATH,"a")
    for read in reads:
        if read.strand == "+":
            sss = read.pos_start_exon+read.pos_end_exon
            sss.sort()
            sss.pop(0)
            if(len(sss+[read.start]) != len(set(sss+[read.start]))):
                continue
            f.write(read.chrom+"\t"+read.name+"\t"+str(read.start))
            for ss in sss[:-1]:
                f.write("\t"+str(ss))
            f.write("\t"+str(read.end))
        else:
            sss = read.pos_start_exon+read.pos_end_exon
            sss.sort()
            sss.pop(len(sss)-1)
            if(len(sss+[read.end]) != len(set(sss+[read.end]))):
                continue
            f.write(read.chrom+"\t"+read.name+"\t"+str(read.end))
            f.write("\t"+str(read.start))
            for ss in sss[1:]:
                f.write("\t"+str(ss))
        f.write("\n")

def visualize_cumulative_sum():
    ssandtss = [ss_mins,tss_mins,ss3end_mins]
    plt.rcParams["font.size"] = 16
    for index in range(3):
        plt.figure()
        dists=[]
        for dis in ssandtss[index]:
            if abs(int(dis)) > args.xlimit:
                dists.append(args.xlimit+2)
                continue
            dists.append(abs(int(dis)))
        cums = stats.cumfreq(dists,numbins=args.xlimit+2)
        plt.xlabel('distance [bp]')
        x = pd.Series(cums.cumcount)
        plt.xlim(0,args.xlimit)
        ax = sns.lineplot(data=x)
        if index==0:
            plt.ylabel('Numper of splice sites')
            plt.savefig("ss_cumulative_plot.png", dpi=500, bbox_inches='tight')
            label="splice site"
        elif index==1:
            plt.ylabel('Numper of TSSs')
            plt.savefig("tss_cumulative_plot.png", dpi=500,bbox_inches='tight')
        else:
            plt.ylabel('Numper of 3-prime ends')
            plt.savefig("ss3end_cumulative_plot.png", dpi=500,bbox_inches='tight')

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="input path for *.bed file", type=str)
    parser.add_argument("ref", help="reference path for gtf or genepred format file", type=str)
    parser.add_argument("--cage", help="optional: path for *.bed cage peaks file", type=str)
    parser.add_argument("--tss_t", help="optional: a threshold for tss correction",type=int, default=40)
    parser.add_argument("--ss_t", help="optional: a threshold for splice site correction", type=int, default=20)
    parser.add_argument("--ss3end_t", help="optional: a threshold for 3'end correction", type=int, default=50)
    parser.add_argument("--xlimit", help="optional: right boundary value for visualization of minimum distance from given elements to annotated values", type=int, default=100)
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = get_args()
    print(args.input)
    ss = input_to_ssdict()
    print("ref:load completed")
    reads = bed.split_per_chrom(bed.read(args.input))
    print("input:load completed")
    ss_mins=[]
    ss3end_mins=[]
    tss_mins=[]
    if hasattr(args,'cage'):
        cage_reads = bed.split_per_chrom(bed.read(args.cage))
        print("cage:load completed")
        reads = bed.compare_per_chrom(reads, cage_reads, find_nearest_CAGEpeak)
        print("cage:correction completed")
    reads = bed.compare_per_chrom(reads, ss, find_nearest_splice_junction)
    print("tss:correction completed")
    bed.execute_per_chrom(reads,write_corrected_ss,OUTPUT_PATH=args.input+".corrected")
    print("corrected file created")
    visualize_cumulative_sum()
