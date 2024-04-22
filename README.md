# exPATT

`exPATT.R` is an R script to detect DNA binding motifs (_e.g._ OCT4) and change the genomic coordinates of each matching position at its center. The script filters genomic regions, recorded in an input BED file, that match a given DNA pattern depicted as a position weight matrix (PWM). These pattern matrices are matched in genomic regions using a minimum threshold given as a percentage of the highest possible score (_e.g._ 80%). Then, genomic regions with a minimum number of hit(s) are filtered and centered around the pattern to generate a new BEDGRAPH file.

## Installation

Hereafter we consider that **R** version 4.3+ is installed on your computer.

Dependencies requieremnts can be installed using the `setup.R` script.
If `Rscript` is in your PATH you just type:

```
> ./setup.R
```

Otherwise you need to source this file from your R IDE:

```
R> source("setup.R")
```

During this installation process various packages will be installed. You can customize the list of packages (especially regarding Genomes annotations from Bioconductor) in `setup.R`.

## Usage

The easiest way to execute `exPATT.R` is directly from your shell command line (given `Rscript` is in your PATH).

A command line helper is shown if no arguments are provided:

```
> exPATT.R
+-----------------------------------------------+
| (c) 2024 - CEA, CNRS, Paris-Saclay University |
+-----------------------------------------------+
$Id: 4309a63ab943ed917080c85d1d33af116620a9d6 $

======================================================================
This script is used to filter a BED file given a set of patterns.

USAGE:
 %prog [-G <genome>] -b <input-bed-file> -p <pattern> [-m <min-score>]
       [-l <min-hits>] [-h <max-hits>] [-w <bed-track-width>]
       -o <output-bed-file>
or
 %prog -L

options (regular mode):
 -G : genome, as a Bioconductor annotation package
      (default: Mmusculus.UCSC.mm9)
 -b : input bed file to be filtered
 -p : the DNA pattern(s), either as a filename or one of the following.
      presets: OCT4,SOX2 (default: OCT4)
 -m : pattern, minimum score (given as percentage of the maximum possible
      score) for counting a match (default: 80)
 -l : lowest/minimum number of hits per BED entries (default:1)
 -h : highest/maximum number of hits per BED entries (default:none)
 -w : width (in bp) of the track in the output file (default: 2000)
 -o : the filename, including the '.bedgraph' extension to use for
      storing results (default: <input-file>_out.bed)

options (genome list mode):
 -L : list all available Bioconductor genomes
======================================================================
```

This program can be executed in two different modes: 
- regular or query mode 
- genomes list, using the `-L` option.

### regular/query mode

This execution mode allows to select regions in an input BED file that match a motif provided as a PWM matrix.

By default, `exPATT` search the _Mmusculus.UCSC.mm9_ genomic regions defined in the input BED file (here `sample.bed`) for the `OCT4` motif with a minimum score of 80% and at least one hit. The output is a BEDGRAPH file of 2kb centered around each hit. The motif is searched on both strands. The execution output is given below: 

```
> exPATT.R -b example/sample.bed
+-----------------------------------------------+
| (c) 2024 - CEA, CNRS, Paris-Saclay University |
+-----------------------------------------------+
$Id: 907226b92e8f6069a12ab8da81c2ef171dca4bef $

- Load Genome: done
- Import BED file: #8432 records
- Retrieve DNA sequences: done
- Search Pattern with min-score=80% (strands '+'/'-'):
  |======================================================================| 100%
  |======================================================================| 100%
- Found 1701 hit(s) on the '+' strand and 1677 on the '-' strand
- Select #hits in the range [1,+Inf]: 1156 on the '+' strand, 1132 on the '-' strand and 545 on both strands
- Filter & Process hit(s) (strands '+'/'-'):
  |======================================================================| 100%
  |======================================================================| 100%
```

`exPATT` includes two default motifs, **OCT4** and **SOX2**. Still, you can use your own motifs using regular tabulated text files to describe the associated PWM matrix. As an example, the PWM matrix of OCT4 is given in the file `example/oct4.pwm`.

### genomes list mode

The `-L` gives a list of all available anotations packages in Bioconductor. Packages already installed in your **R** environment are indicated with a star character (*). For example, the genome annotations _Mmusculus.UCSC.mm8_ is installed as shown below.

```
> ./exPATT.R -L
+-----------------------------------------------+
| (c) 2024 - CEA, CNRS, Paris-Saclay University |
+-----------------------------------------------+
$Id: 907226b92e8f6069a12ab8da81c2ef171dca4bef $


Genome                                   :                       OS       Prov.                Ver. Inst.
Alyrata.JGI.v1                           :                  Alyrata         JGI                  v1
Amellifera.BeeBase.assembly4             :               Amellifera     BeeBase           assembly4
Amellifera.NCBI.AmelHAv3.1               :               Amellifera        NCBI          AmelHAv3.1
...
Mmusculus.UCSC.mm39                      :                Mmusculus        UCSC                mm39
Mmusculus.UCSC.mm8                       :                Mmusculus        UCSC                 mm8 *
Mmusculus.UCSC.mm8.masked                :                Mmusculus        UCSC                 mm8
...
Vvinifera.URGI.IGGP12Xv0                 :                Vvinifera        URGI           IGGP12Xv0
Vvinifera.URGI.IGGP12Xv2                 :                Vvinifera        URGI           IGGP12Xv2
Vvinifera.URGI.IGGP8X                    :                Vvinifera        URGI              IGGP8X
```

To install new annotations package you must do it from an **R** session. For example, if you want to install the _Mmusculus.UCSC.mm39_ version you must use the following command on the R console.

```
R> BiocManager::install("BSgenome.Mmusculus.UCSC.mm39")
```

Notice the **BSgenome** prefix that is not shown on the genomes list to optimize the display.

## Reference

