# long_read_isoform_visualizer
Simple sashimi-plot like isoform visualizer for long-read sequencing.



## Features
- Plots isoforms from full-length transcripts genereted by long-read sequencer such as ONT MinION and PacBio.
- Plots RNAseq reads densities along exons and splice junctions.
- Plots isoforms transcribed from given Transcript Start Sites (TSSs)
    - Allows multiple TSSs on same figure and can compare isoforms from different TSSs.
- Plots all isoforms in the given regions.

## Requirement

- Python3 (3.7.6) 

## Reference files
We checked performance using following reference files.

- FANTOM5 CAGE peaks
    - http://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/CAGE_peaks/
- GENCODE annotation files (We used ver 30)
    - https://www.gencodegenes.org/human/
- RefSeq
    - http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/

## Installation

```
git clone https://github.com/isihya/long_read_isoform_visualizer.git
```

## Usage
First, please conduct correction by reference files.

```
usage: correction.py [-h] [--cage CAGE] [--tss_t TSS_T] [--ss_t SS_T]
                     [--ss3end_t SS3END_T] [--xlimit XLIMIT]
                     input ref

positional arguments:
  input                input path for *.bed file
  ref                  reference path for gtf or genepred format file

optional arguments:
  -h, --help           show this help message and exit
  --cage CAGE          optional: path for *.bed cage peaks file
  --tss_t TSS_T        optional: a threshold for tss correction
  --ss_t SS_T          optional: a threshold for splice site correction
  --ss3end_t SS3END_T  optional: a threshold for 3'end correction
  --xlimit XLIMIT      optional: right boundary value for visualization of
                       minimum distance from given elements to annotated
                       values
```

Then, you can plot isoforms.

```
usage: plot.py [-h] [--tss TSS] [--cellname CELLNAME] [--xliml XLIML]
               [--xlimr XLIMR]
               input chrm

positional arguments:
  input                path for *.corrected file
  chrm                 chromosome. set 'chrN' (N=1,...,22,X,Y)

optional arguments:
  -h, --help           show this help message and exit
  --tss TSS            optional: coordinate of tss which you want to extract
                       by
  --cellname CELLNAME  optional: you can add cellname on figure
  --xliml XLIML        optional: left boundary value for visualization
  --xlimr XLIMR        optional: right blundary value for visualization
```
