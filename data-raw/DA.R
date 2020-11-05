#!/bin/r

if(sub('.*\\/','',getwd())=='data-raw'){
	 setwd('..')
}

library(DBI)
library(DESeq2)
library(dirfns)

source("R/DESeqFns.R")
source("R/sqlfns.R")

con <- dbConnect(RSQLite::SQLite(),'data/atacCiona.db')

expDesign <- dbReadTable(con,'ataclib',row.names='lib')
# expDesign <- read.delim('data-raw/dat/meta/ataclib.txt',row.names=1)
# expDesign <- expDesign[!as.logical(expDesign$omit),]

# read counts from htseq
counts <- dbReadTable(con,'atacreads',row.names='PeakID')
# ATAC-seq counts
# counts <- lrtab(
#   'data-raw/dat/counts/',
#   read.table
# )
# dat <- counts[-grep('rnaseq',names(counts))]
# dat <- do.call(cbind,sapply(dat,'[',4))
# colnames(dat) <- make.names(sub('_counts.V4','',colnames(dat)))
# row.names(dat) <- peak.id.coord[do.call(function(...)paste(...,sep = ':'),counts[[1]][,1:3]),]
# dat <- dat[!is.na(row.names(dat)),]
# 
# dat <- as.data.frame(dat)
# dat <- dat[,ataclib[,1]]

filt <- counts[,row.names(expDesign)]

# remove replicates with fewer than 500000 reads
filt <- filt[,apply(filt,2,sum)>500000]
# reorder expDesign by column order of filt
expDesign <- expDesign[names(filt),]

# tissue-dependent model
design.gfp <- subset(expDesign,condition=='control'&tissue!='Ef1a'&time!='6hpf'&!KO)
dat.gfp <- filt[,row.names(design.gfp)]
gfp <- get.dds(
  filt,design.gfp,~0+time*tissue+method,
  path='gfp',
  list(c('tissue','B7.5','mesenchyme')),
  test="LRT",reduced=~0+time+tissue+method,cooksCutoff=F,alpha=.05
)

# subset peaks by B7.5 lineage
b75expDesign <- subset(expDesign,tissue=='B7.5')
b75 <- filt[,row.names(b75expDesign)]
filt <- filt[row.names(b75),]

# 10hpf model
design10 <- subset(b75expDesign,time=='10hpf')
mesp <- get.dds(
  filt[,row.names(design10)],design10,~0+condition+KO,
  path='cond10',
  comparisons=list(
   c('condition','FoxF_KO',"control"),
    c('condition','mesp_MekMut',"mesp_dnFGFR"),
    c('condition','mesp_MekMut',"control"),
    c('condition','mesp_dnFGFR',"control")
  ),
  test='LRT',reduced=~0+KO,cooksCutoff=F,alpha=.05
)

# 15-20hpf model
designHandr <- subset(b75expDesign,time%in%c('15hpf','18hpf','20hpf'))
handr <- get.dds(
  filt[,row.names(designHandr)],designHandr,~0+condition+time+method,
  path='handr',
  comparisons=list(
    c('condition','handr_MekMut',"handr_dnFGFR"),
    c('condition','handr_MekMut',"control"),
    c('condition','handr_dnFGFR',"control")
  ),
  test='LRT',reduced=~0+time+method,cooksCutoff=F,alpha=.05
)

# time-dependent model
designLacz <- subset(b75expDesign,condition=='control'&!KO)
lacz <- get.dds(
  filt[,row.names(designLacz)],designLacz,~0+time+method,
  path='time',
  # pairwise comparisons for LRT
  comparisons=list(
    c('time','6hpf',"10hpf"),
    c('time','10hpf',"15hpf"),
    c('time','10hpf',"18hpf")
  ),
  test='LRT',reduced=~0+method,cooksCutoff=F,alpha=.05
)

peakdat <- lapply(
  list(mesp,handr,lacz,gfp),
  function(x) counts(x$dds,normalized=T)
)

mapply(
  dbWritePeaks,
  name=c('mespct','handrct','timect','tissuect'),
  peakdat,
  MoreArgs = list(conn=con,row.names="PeakID",overwrite=T)
)

atac <- sapply(
  Reduce(
    function(x,y) append(x, y$res),
    list(mesp,handr,lacz,gfp),
    init = NULL
  ),
  function(z) cbind(
    PeakID=row.names(z),
    as.data.frame(z),
    stringsAsFactors=F
  ),
  simplify = F
)

dbWritePeaks(con,'atacseq',melt.rename(atac,'comparison'))
