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

## Installation

```
git clone https://github.com/isihya/long_read_isoform_visualizer.git
```

## Usage
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
