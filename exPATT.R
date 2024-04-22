#!Rscript --vanilla
#
# exPATT
# Copyright (C) 2024  CEA, CNRS, Paris-Saclay University
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

VERSION <- "$Id$"

# Constants
OCT4 <- structure(c(0.007407407, 0.103703704, 0.4, 0.011111111, 0.274074074, 
                    0.144444444, 0.022222222, 0, 0.011111111, 0.848148148, 
                    0.014814815, 0.011111111, 0.051851852, 0.3, 0, 0.055555556,
                    0.022222222, 0.007407407, 0.055555556, 0.07037037, 
                    0.137037037, 0.17037037, 0, 0.025925926, 0.02962963, 
                    0.044444444, 0.022222222, 0.011111111, 0.851851852, 
                    0.025925926, 0, 0.007407407, 0.007407407, 0.992592593, 
                    0, 0, 0.022222222, 0, 0.007407407, 0.896296296, 0.003703704, 
                    0.07037037, 0.314814815, 0.037037037, 0.422222222, 
                    0.092592593, 0.003703704, 0.014814815, 0.37037037, 0, 
                    0.022222222, 0.1, 0.807407407, 0.003703704, 0.018518519, 
                    0.003703704, 0.040740741, 0.233333333, 0.011111111, 
                    0.044444444), 
                  .Dim = c(4L, 15L), 
                  .Dimnames = list(c("A", "T", "C", "G"), NULL)
)
OCT4 <- OCT4[c("A","C","G","T"),]

SOX2 <- structure(c(0.018087855, 0.196382429, 0.015503876, 0.087855297, 
                    0.167958656, 0.005167959, 0.005167959, 0.160206718, 
                    0.043927649, 0.002583979, 0, 0.842377261, 0.07751938, 0, 
                    0.018087855, 0.625322997, 0.025839793, 0.056847545, 
                    0.188630491, 0.018087855, 0.124031008, 0.002583979, 
                    0.018087855, 0.30749354, 0.015503876, 0.018087855, 
                    0.015503876, 0.03875969, 0.813953488, 0.028423773, 
                    0.007751938, 0.007751938, 0.007751938, 0.007751938, 0, 
                    0.984496124, 0.036175711, 0.007751938, 0.813953488, 
                    0.005167959, 0.007751938, 0.248062016, 0.051679587, 
                    0.049095607, 0.470284238, 0.005167959, 0.007751938, 
                    0.124031008, 0.235142119, 0.020671835, 0.085271318, 
                    0.012919897, 0.599483204, 0.012919897, 0.033591731, 
                    0.046511628, 0.031007752, 0.028423773, 0.033591731, 
                    0.26873385), 
                  .Dim = c(4L, 15L), 
                  .Dimnames = list(c("A", "C", "G", "T"),NULL))

SOX2 <- SOX2[c("A","C","G","T"),]


PATTERNS <- c("OCT4","SOX2")

cat("+-----------------------------------------------+\n")
cat("| (c) 2024 - CEA, CNRS, Paris-Saclay University |\n")
cat("+-----------------------------------------------+\n")
cat(VERSION,"\n\n")

options(stringsAsFactors=FALSE)
suppressPackageStartupMessages(library(tools))

## FUNCTIONS
## ----------------------------------------

init_libs <- function() {
  options(warn=-1)
  library(parallel)
  suppressPackageStartupMessages(library(ngs,quietly=TRUE,verbose=FALSE,warn.conflicts=FALSE))
  suppressPackageStartupMessages(library(Biostrings,quietly=TRUE,verbose=FALSE,warn.conflicts=FALSE))
  suppressPackageStartupMessages(library(BSgenome,quietly=TRUE,verbose=FALSE,warn.conflicts=FALSE))
  suppressPackageStartupMessages(library(rtracklayer,quietly=TRUE,verbose=FALSE,warn.conflicts=FALSE))
  suppressPackageStartupMessages(library(GenomicRanges,quietly=TRUE,verbose=FALSE,warn.conflicts=FALSE))
  options(warn=0)
}

print_help <- function() {
  cat("======================================================================\n")
  cat("This script is used to filter a BED file given a set of patterns.\n\n")
  cat("USAGE:\n")
  cat(" %prog [-G <genome>] -b <input-bed-file> -p <pattern> [-m <min-score>]\n")
  cat("       [-l <min-hits>] [-h <max-hits>] [-w <bed-track-width>]\n")
  cat("       -o <output-bed-file>\n")
  cat("or\n")
  cat(" %prog -L\n\n")
  cat("options (regular mode):\n")
  cat(" -G : genome, as a Bioconductor annotation package\n")
  cat("      (default: Mmusculus.UCSC.mm9)\n")
  cat(" -b : input bed file to be filtered\n")
  cat(" -p : the DNA pattern(s), either as a filename or one of the following.\n")
  cat("      presets: OCT4,SOX2 (default: OCT4)\n")
  cat(" -m : pattern, minimum score (given as percentage of the maximum possible\n") 
  cat("      score) for counting a match (default: 80)\n")
  cat(" -l : lowest/minimum number of hits per BED entries (default:1)\n")
  cat(" -h : highest/maximum number of hits per BED entries (default:none)\n")
  cat(" -w : width (in bp) of the track in the output file (default: 2000)\n")
  cat(" -o : the filename, including the '.bedgraph' extension to use for\n")
  cat("      storing results (default: <input-file>_out.bed)\n\n")
  cat("options (genome list mode):\n")
  cat(" -L : list all available Bioconductor genomes\n")
  cat("======================================================================\n\n")
}

parse_cmdline <- function(args) {
  genome <- "Mmusculus.UCSC.mm9"
  pattern <- "OCT4"
  pmin_score <- "80%"
  min_hits <- 1
  max_hits <- NA
  bed_width <- 2000
  bed_file <- NA
  out_file <- NA
  list_mode <- FALSE
  state <- "X"
  for(i in 1:length(args)) {
    opt <- args[[i]]
    if(substring(opt,1,1) == "-") {
      m  <- substring(opt,2,2)
      state <- switch(m,g="G",b="b",p="p",m="m",o="o",l="l",h="h",w="w",L="L","X")
      if(state=="L") return(list(list_mode=TRUE))
    }
    else {
      if(state == "G")
        genome <- opt
      else if(state=="b")
        bedfile <- opt
      else if(state=="p")
        pattern <- opt
      else if(state=="m")
        pmin_score <- paste(opt,"%",sep="")
      else if(state=="l")
        min_hits <- as.integer(opt)
      else if(state=="h")
        max_hits <- as.integer(opt)
      else if(state=="w")
        bed_width <- as.integer(opt)
      else if(state=="o")
        out_file <- opt      
      else
        stop("Invalid command line..")
    }
  }
  
  # check bed-width
  if(is.na(bed_width)) {
    printHelp()
    cat("\nERROR: Invalid bed track width...\n")
    quit()  
  }
  
  # check bed-file
  if(is.na(bedfile)) {
    printHelp()
    quit()
  }
  
  # check min-hits
  if(min_hits < 1) {
    printHelp()
    cat("\nERROR: <min-hits> must be greater or equal to 1\n")
    quit()
  }
  
  # check max-hits
  if(!is.na(max_hits) && max_hits < min_hits) {
    printHelp()
    cat("\nERROR: <max-hits> must be greater or equal to <min-hits>\n")
    quit()
  }
  
  # set pattern
  if(pattern %in% PATTERNS) {
    pattern <- get(pattern)
  }
  else {
    if(!file.exists(pattern)) stop("Pattern file doesn't exists...")
    pattern <- as.matrix(read.delim2(pattern,
                                     header = FALSE,
                                     row.names = 1,
                                     dec = "."))
    if(!all(row.names(pattern) == c("A","C","G","T")))
      stop("Invalid PWM Pattern file...", row.names(pattern))
  }
  
  # set output file
  if(is.na(out_file))
    out_file <- gsub(".bed","_out.bed",bedfile)
  
  # eop
  return(list(genome=genome,pattern=pattern,pmin_score=pmin_score,
              min_hits=min_hits, max_hits=max_hits, bed_width=bed_width,
              bedfile=bedfile, output=out_file))
}

MAIN <- function(args) {
  
  # init.
  pat.delta <- ncol(args$pattern) %/% 2
  
  # load genome
  genome <- paste("BSgenome", args$genome, sep=".")
  cat("- Load Genome: ")
  LoadGenome(genome)
  cat("done\n")
  
  # load the BED source file
  cat("- Import BED file: ")
  bed.data <- import(args$bedfile)
  cat("#",length(bed.data)," records\n",sep="")
  
  # @TODO: la normalisation doit avoir lieu ICI
  
  # get sequences
  cat("- Retrieve DNA sequences: ")
  dna.seq <- getSeq(get(genome),bed.data)
  names(dna.seq) <- as.character(seqnames(bed.data))
  assign("dna.seq",dna.seq,envir=globalenv())
  cat("done\n")
  
  # match patterns
  cat("- Search Pattern with min-score=",args$pmin_score,
      " (strands '+'/'-'):\n",sep="")
  nb.seqs <- length(dna.seq)
  cpt <- 1
  pb <- txtProgressBar(min = 0, max = nb.seqs, initial = cpt, style = 3)
  hits <- sapply(dna.seq, function(s) {
    h <- matchPWM(args$pattern, s, min.score = args$pmin_score)
    cpt <<- cpt + 1
    setTxtProgressBar(pb,cpt)
    return(h)
  })
  close(pb)
  
  cpt <- 1
  pb <- txtProgressBar(min = 0, max = nb.seqs, initial = cpt, style = 3)
  hits.rev <- sapply(dna.seq,function(s) {
    h <- matchPWM(reverseComplement(args$pattern),s)
    cpt <<- cpt + 1
    setTxtProgressBar(pb,cpt)
    return(h)
  })
  close(pb)
  
  nhits <- sapply(hits,length)
  nhits.rev <- sapply(hits.rev,length)
  fhits <-  nhits > 0
  fhits.rev <- nhits.rev > 0
  cat("- Found ",sum(fhits)," hit(s) on the '+' strand and ", sum(fhits.rev), 
      " on the '-' strand\n",sep="")
  
  # apply min/max hit bounds
  minH <- args$min_hits
  maxH <- args$max_hits
  if(is.na(maxH)) maxH <- max(nhits+nhits.rev) + 1
  flag <- (nhits + nhits.rev) > maxH
  nhits[flag] <- 0
  nhits.rev[flag] <- 0
  nhits[nhits < minH | nhits > maxH] <- 0
  nhits.rev[nhits.rev < minH | nhits.rev > maxH] <- 0
  fhits <-  nhits > 0
  fhits.rev <- nhits.rev > 0
  cat("- Select #hits in the range [",args$min_hits,",",
      ifelse(is.na(args$max_hits),"+Inf",args$max_hits),"]: ",
      sum(fhits & !(fhits & fhits.rev))," on the '+' strand, ", 
      sum(fhits.rev & !(fhits & fhits.rev)), 
      " on the '-' strand and ",sum(fhits & fhits.rev),
      " on both strands\n",sep="")
  
  # '+' strand - filter & hits
  bed.plus <- bed.data[fhits]
  hits.plus <- hits[fhits]
  cat("- Filter & Process hit(s) (strands '+'/'-'):\n")
  
  # -- refactor as a function
  cpt <- 1
  pb <- txtProgressBar(min=0, max=length(bed.plus),initial=cpt, style=3)
  GR.plus <- sapply(1:length(bed.plus), function(idx) {
    setTxtProgressBar(pb,cpt)
    bed <- bed.plus[idx]
    meta <- elementMetadata(bed)
    hits <- hits.plus[[idx]]
    chr <-rep(seqnames(bed), length.out = length(hits))
    strand <- rep(strand(bed),length.out = length(hits))
    ir <- sapply(1:length(hits), function(i) {
      h <- hits[i]
      c <- start(ranges(bed)) + start(h) + pat.delta - 1 
      start <- c - args$bed_width
      end <- c + args$bed_width
      return(c(start,end))
    })
    gr <- GRanges(seqnames=chr, 
                  ranges=IRanges(start=ir[1,],end=ir[2,]),
                  strand=strand)
    elementMetadata(gr) <- meta
    cpt <<- cpt + 1
    return(gr)
  })
  GR.plus <- do.call(c,GR.plus)
  close(pb)
  
  # '-' strand - filter & hits
  bed.minus <- bed.data[fhits.rev]
  hits.minus <- hits.rev[fhits.rev]
  
  cpt <- 1
  pb <- txtProgressBar(min=0, max=length(bed.minus),initial=cpt, style=3)
  GR.minus <- sapply(1:length(bed.minus), function(idx) {
    setTxtProgressBar(pb,cpt)
    bed <- bed.minus[idx]
    meta <- elementMetadata(bed)
    hits <- hits.minus[[idx]]
    chr <-rep(seqnames(bed), length.out = length(hits))
    strand <- rep(invertStrand(strand(bed)),length.out = length(hits))
    ir <- sapply(1:length(hits), function(i) {
      h <- hits[i]
      c <- start(ranges(bed)) + start(h) + pat.delta - 1 
      start <- c - args$bed_width
      end <- c + args$bed_width
      return(c(start,end))
    })
    gr <- GRanges(seqnames=chr, 
                  ranges=IRanges(start=ir[1,],end=ir[2,]),
                  strand=strand)
    elementMetadata(gr) <- meta
    cpt <<- cpt + 1
    return(gr)
  })
  GR.minus <- do.call(c,GR.minus)
  close(pb)
  
  # export results (rtracklayer)
  export(c(GR.plus,GR.minus),args$output, format="bed")
  
  # eop
  return(invisible())
}

ListGenomes <- function() {
  
  # get genomes
  genomes <- available.genomes(TRUE)
  installed <- installed.genomes()
  
  # add installed status
  genomes$installed <- sapply(genomes$pkgname, 
                              function(x) {
                                if(x %in% installed)
                                  return("*")
                                else
                                  return ("")
                              })
  
  # filter genomes names
  genomes$pkgname <- sapply(genomes$pkgname,function(x) {gsub("BSgenome\\.","",x)})
  
  # build formatting string
  len_pkgname <- max(nchar(genomes$pkgname))
  len_os <- max(nchar(as.character(genomes$organism)))
  len_provider <- max(nchar(as.character(genomes$provider)))
  len_version <- max(nchar(as.character(genomes$genome)))
  format_string <- sprintf("%%-%ds: %%%ds %%%ds %%%ds %%s\n",
                           as.integer(len_pkgname),
			   as.integer(len_os),
			   as.integer(len_provider),
			   as.integer(len_version))
  
  # display results
  cat("\n")
  cat("\033[1m",
      sprintf(format_string, 
              "Genome", "OS", "Prov.", "Ver.", "Inst."),
      "\033[0m", sep="")
  tmp <- apply(genomes, 1,
                function(x,f) {
                  cat(sprintf(f,x["pkgname"],trimws(x["organism"]),x["provider"],
                                x["genome"],x["installed"]))
                }, format_string)
  cat("\n")
}

LoadGenome <- function(genome) {

  # init.
  avail.genomes <- installed.genomes()
  genome.pkg <- paste(genome,sep=".")
  
  # check
  if(!genome.pkg %in% avail.genomes) {
    cat("<not-found>\n")
    cat("- Install new genome package: ",genome.pkg,"\n\n")
    BiocManager::install(genome.pkg)
  }
  
  # load
  do.call("library",list(package=genome.pkg))
  
  # eop
  return(invisible(TRUE))
}

## MAIN
## ----------------------------------------

if(!interactive()) {
  
  args <- commandArgs(trailingOnly=TRUE)
  
  if(length(args) == 0) {
    print_help()  
    quit()
  }
  
  init_libs()
  
  args <- parse_cmdline(args)
  
  if("list_mode" %in% names(args)) {
    ListGenomes()
  } else {
    
    if(length(args$pattern) == 0) {
      print_help()
      stop("The 'pattern' argument is required\n")
    }
    
    MAIN(args)
  }
}