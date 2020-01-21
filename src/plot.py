

import bed_utils as bed
from collections import Counter
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from pylab import rcParams
from scipy import signal, interpolate
import numpy as np
import inline
import argparse

# 最も多いcoverageの高さ=グラフの高さのmaxとし、領域ごとに相対的な高さを求める
# 領域ごとにリードを長方形で表現
# 3次スプライン曲線で、長方形をつなぐ。
# スプライン曲線の高さは最も多いcoverageの高さ/2とする
# スプライン曲線の頂点に白抜きの正方形を上から表示、その部分に数を表現


# 各splice site間のつながりの数を求める
# 各splice site間のcoverageの深さを求める


def input(INPUT_PATH):
    f = open(INPUT_PATH)
    reads = f.readlines()
    return reads

def preprocess_from_TSS(reads,TARGET_TSS):
    field = Counter()
    for read in reads:
        read = read.rstrip("\n").split("\t")
        if(read[0]!=args.chrm):
            continue
        read = read[1:]
        if(int(read[1]) == TARGET_TSS):
            #if len(read) == 3:
            #    continue
            readname = read[0]
            ss = read[1:]
            ss.sort()
            #print(ss)
            flag = 0
            for i in range(len(ss)-1):
                #print([ss[i],ss[i+1],str(flag%2)])
                field.update([",".join([ss[i],ss[i+1],str(flag%2)])])
                flag += 1
    return field


def preprocess_from_Region(reads,xliml,xlimr):
    field = Counter()
    for read in reads:
        read = read.rstrip("\n").split("\t")
        if(read[0]!=args.chrm):
            continue
        read = read[1:]
        ss = read[1:]
        ss.sort()
        if int(ss[0]) > xlimr or int(ss[-1]) < xliml :
            continue
        flag = 0
        for i in range(len(ss)-1):
            #print([ss[i],ss[i+1],str(flag%2)])
            field.update([",".join([ss[i],ss[i+1],str(flag%2)])])
            flag += 1
    return field



def make_name_to_read_end_dict(INPUT_PATH):
    reads = bed.read(INPUT_PATH)
    nametoend = {}
    for read in reads:
        nametoend[read.name] = read.end
    return nametoend

def calculate_height(field,TARGET_TSS,color,offset,ax):
    # flag = 0 means exon
    sslist = []
    for sstuple in field.keys():
        sstuple = sstuple.split(",")
        sstuple = list(map(int,sstuple))
        if sstuple[0] not in sslist:
            sslist.append(sstuple[0])
        if sstuple[1] not in sslist:
            sslist.append(sstuple[1])
    sslist.sort()
    hist = [0]*len(sslist)
    flag = 0
    for ss1 in range(len(sslist)):
        for ss2 in range(len(sslist)):
            if ss1 >= ss2:
                continue
            if(flag == 0):
                key = ",".join([str(sslist[ss1]),str(sslist[ss2]),"0"])
                hist[ss1]+=field[key]
                hist[ss2]-=field[key]
    for i in range(len(hist)-1):
        hist[i+1]+=hist[i]
    if(offset == [-1]):
        offset = [0]*len(sslist)
    if(sslist==[]):
        return ax
    ymax=0
    for i in range(len(sslist)-1):
        for j in range(len(sslist)-1):
            if(j<i):
                continue
            if args.xliml!=None and args.xlimr!=None:
                if args.xliml>sslist[i] or args.xlimr<sslist[j]:
                    continue
            width=sslist[j]-sslist[i]
            ## splice sitesをつなぐ線
            x = np.array([sslist[i],sslist[i]+(sslist[j]-sslist[i])/2,sslist[j]])
            #y = np.array([hist[i],max(hist[i],hist[i+1])/2,hist[i+1]])
            key = ",".join([str(sslist[i]),str(sslist[j]),"1"])
            if(field[key]< 1):
                continue
            height = field[key]
            if(field[key]<10):
                height = field[key]*13
            elif(field[key]<100):
                height = field[key]*4
            ymax=max(height,ymax)
            y = np.array([0,height,0])
            f = interpolate.interp1d(x, y, kind="quadratic")
            xnew = np.linspace(sslist[i],sslist[j],51)
            ax.plot(xnew, f(xnew), '-',lw=0.5, color="black", zorder=20)
            plt.text(sslist[i]+(sslist[j]-sslist[i])/2, height, str(field[key]), zorder=30 , fontsize=10)
    for i in range(len(sslist)-1):
        width=sslist[i+1]-sslist[i]
        ymax=max(ymax,offset[i]+hist[i])
        r = patches.Rectangle(xy=(int(sslist[i]), 0), width=width, height=offset[i]+hist[i], lw=0.5, ec="black", color=color, zorder=1)
        offset[i]+=hist[i]
        ax.add_patch(r)
    ax.ticklabel_format(style='sci',axis='x',scilimits=(0,0))
    #plt.title('All isoforms transcribed from chr'+str(chrm)+": "+str(TARGET_TSS)+", cell line: "+args.cellname)
    #cell-lineは図のcaptionにつけることにする
    if TARGET_TSS != "":
        plt.title('All isoforms transcribed from '+str(args.chrm)+", "+str(TARGET_TSS))
    else:
        plt.title('All isoforms from '+str(args.chrm)+": "+str(args.xliml)+" to "+str(args.xlimr))
    plt.ylabel('Count of reads')
    if args.xliml!=None and args.xlimr!=None :
        plt.xlim([args.xliml,args.xlimr])
    elif args.xliml==None and args.xlimr==None:
        pass
    else:
        raise "please set xliml and xlimr together"

    plt.ylim([0,ymax*1.1])

    #plt.savefig(args.cellname+",TSS:"+str(TARGET_TSS)+".png",dpi=500)
    return ax

# isoformごとに色を分けたい
def preprocess_per_isoform(reads,TARGET_TSS):
    isoforms = {}
    for read in reads:
        read = read.rstrip("\n").split("\t")[1:]
        if(int(read[1]) == TARGET_TSS and len(read)>3):
            if(int(read[1])>int(read[2])):
                if(",".join(read[2:]) not in isoforms):
                    isoforms[ ",".join(read[2:])]=[]
                else:
                    isoforms[ ",".join(read[2:])].append(read)
            else:
                if(",".join(read[1:-1]) not in isoforms):
                    isoforms[ ",".join(read[1:-1])]=[]
                else:
                    isoforms[ ",".join(read[1:-1])].append(read)
    return isoforms

def make_color():
    palette = ["r", "g", "b", "c", "m", "y","lightskyblue","darkviolet","navy","chocolate"]
    return palette

def print_all():
    offset = [-1]
    if args.tss!=None:
        for i in range(len(TARGET_TSSs)):
            print(str(i+1)+" of "+ str(len(TARGET_TSSs)) + " TSS processing")
            TARGET_TSS = int(TARGET_TSSs[i])
            TARGET_TSS = int(TARGET_TSSs[i])
            field = preprocess_from_TSS(reads,TARGET_TSS)
            ax = plt.subplot(len(TARGET_TSSs),1,i+1)
            ax = calculate_height(field,TARGET_TSS,"skyblue",offset,ax)
            if(i==len(TARGET_TSSs)-1):
                ax.set_xlabel('Coordinates on hg38')
            else:
                ax.tick_params(labelbottom=False)
    else:
        field = preprocess_from_Region(reads,args.xliml,args.xlimr)
        fig = plt.figure(figsize=(8.0, 4.0))
        ax = fig.add_subplot(111)
        ax = calculate_height(field,"","skyblue",offset,ax)
        ax.set_xlabel('Coordinates on hg38')
    plt.subplots_adjust(hspace=0.5)
    #plt.tight_layout()
    if len(args.cellname) > 0:
        cellnameoffset=args.cellname+","
    else:
        cellnameoffset=""
    if args.tss!=None:
        plt.savefig("longread_isoforms_"+cellnameoffset+"TSS:"+str(TARGET_TSS)+".png",dpi=500, bbox_inches="tight")
    else:
        plt.savefig("longread_isoforms_"+cellnameoffset+"Region:"+str(args.xliml)+"-"+str(args.xlimr)+".png",dpi=500, bbox_inches="tight")



def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="please set me", type=str)
    parser.add_argument("chrm", help="please set me", type=str)
    parser.add_argument("--tss", help="option", type=str)
    parser.add_argument("--cellname", help="option", type=str, default="")
    parser.add_argument("--xliml", help="option", type=int)
    parser.add_argument("--xlimr", help="option", type=int)
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = get_args()
    print(args.input)
    if args.tss==None and (args.xliml==None and args.xlimr==None):
        raise "please set tss or {xlim,xlimr}"
    plt.rcParams["font.size"] = 13
    if args.tss is not None:
        TARGET_TSSs=args.tss.split(",")
    reads = input(args.input)
    palette = make_color()
    print_all()
