#!/bin/r

if(sub('.*\\/','',getwd())=='data-raw'){
	 setwd('..')
}

# initialize the database used for all subsequent analysis
library(DBI)
library(GenomicFeatures)
library(rtracklayer)
library(BSgenome.Crobusta.HT.KY)
library(peakToGene)
library(dirfns)

source('R/sqlfns.R')

# initialize database
con <- dbConnect(RSQLite::SQLite(),'data/atacCiona.db')

# get KH-to-KY table
kh.ky <- read.table('data-raw/KH2012_KY2019.txt',stringsAsFactors = F)
names(kh.ky) <- c('KHID',"KYID")
kh.ky$KHID <- paste0("KH2013:",kh.ky$KHID)
kh.ky$KYID <- paste0("KY2019:",kh.ky$KYID)
dbWriteTable(con, "KH_KY", kh.ky)

# read gene data
genedat <- lrtab("data-raw/dat/genes/",read.delim,'txt',quote='',stringsAsFactors=F)
genedat$ATM_genes_from_ANISEED <- genedat$ATM_genes_from_ANISEED[!duplicated(genedat$ATM_genes_from_ANISEED),]
genedat <- lapply(genedat, function(x) setNames(x,sub("GeneID", "KHID", names(x))))

htky <- getFeatures('data-raw/HT.Gene.gff3')

# TSS-seq
tsc <- read.csv('data-raw/dat/peaks/tsc.csv')
tsc$KHID <- paste0("KH2013:",tsc$GeneID)
dbWriteKey(con, 'tsc', tsc, foreign = 'KH_KY(KHID)')

tsc <- GRanges(
  Rle(tsc$Chr),
  IRanges(tsc$Start,tsc$End),
  Rle(tsc$Strand),
  KHID=tsc$KHID,
  Rep.TSS=tsc$Rep.TSS,
  feature=tsc$Location
)

# annotate peaks
peakome <- import('data-raw/accessomeKY.bed')
peakome <- peakome[width(peakome)>50]
names(peakome) <- peakome$name

# promoterAnn <- findOverlaps(peakome,htky$promoter)

genomefeat <- append(
  #lapply(htky$features,unlist),
  htky, list(TSC=tsc,genome=GRanges(seqinfo(Crobusta)))
)

overlaps <- lapply(htky, getOverlaps, peakome)
overlaps <- lapply(overlaps, setNames, c("PeakID", "KYID"))
# genomefeat <- reduce(do.call(GRangesList,genomefeat))
# genomefeat <- sapply(genomefeat,trim)

# peakGeneAnnotation <- annotatePeaks(peakome,htky$genebody,features = genomefeat)
# peakGeneAnnotation$promoterAnn <- promoterAnn

peakfeat <- sapply(overlaps, function(x) peakome$name%in%x[,1])
row.names(peakfeat) <- peakome$name
dbWriteKey(con, 'peakfeature', peakfeat*1 , primary = 'PeakID', row.names = "PeakID")

# gene.peak <- reshape2::melt(peakGeneAnnotation$peakToGene)
gene.peak <- do.call(rbind, overlaps)
gene.peak <- gene.peak[!duplicated(gene.peak),]
# names(gene.peak) <- c("GeneID","PeakID")
dbWriteKey(con,'geneToPeak',gene.peak,foreign = c("KH_KY(KYID)","peakfeature(PeakID)"))

# genebody <- data.frame(
#   GeneID=names(peakGeneAnnotation$genes),
#   chr=seqnames(peakGeneAnnotation$genes),
#   start=start(peakGeneAnnotation$genes)+10000,
#   end=end(peakGeneAnnotation$genes)-10000,
#   stringsAsFactors = F
# )
# add.genes <- genebody$GeneID[!genebody$GeneID%in%genedat$gene_name$GeneID]
# genedat$gene_name[add.genes,] <- rep(add.genes,2)
# genedat$gene_name <- merge(genebody,genedat$gene_name)
# dbWriteKey(con, "gene_name", genedat$gene_name, primary = "GeneID")

mapply(
  dbWriteKey,
  name=names(genedat),
  genedat,
  MoreArgs = list(conn=con, foreign = "KH_KY(KHID)")
)

# scRNAseq
scrna <- lrtab('data-raw/dat/scrna/',read.csv,stringsAsFactors=F)
scrna <- lapply(scrna, function(x) setNames(x,sub("GeneID", "KHID", names(x))))

scrna$TVCP <- as.data.frame(rbind(as.matrix(scrna$TVCP),cbind(
  c("KH2013:KH.L71.1","KH2013:KH.C10.370","KH2013:KH.C1.404","KH2013:KH.C1.279"),
  dbReadTable(con,'gene_name',row.names="KHID")[c(
    "KH2013:KH.L71.1","KH2013:KH.C10.370","KH2013:KH.C1.404","KH2013:KH.C1.279"
  ),"UniqueNAME"],"TVC-specific","Primed"
)))

dbWriteKey(con, 'scrna', melt.rename(scrna,'geneset'), foreign='KH_KY(KHID)')

ebfdat <- lrtab('data-raw/dat/ebfdat/',read.csv,stringsAsFactors=F)
ebfdat <- lapply(ebfdat, function(x) setNames(x,sub("GeneID", "KHID", names(x))))
dbWriteKey(con, 'ebfdat', melt.rename(ebfdat,'geneset'), foreign="KH_KY(KHID)")

genedat <- lrtab("data-raw/dat/genes/",read.csv,'csv',stringsAsFactors=F)
genedat <- lapply(genedat, function(x) setNames(x,sub("GeneID", "KHID", names(x))))

names(genedat$mesenchyme20hpf)[1] <- "UniqueNAME"
dbWriteKey(con,'mesenchyme20hpf',genedat$mesenchyme20hpf,foreign = "gene_name(UniqueNAME)")

genedat$mesenchyme20hpf <- do.call(data.frame,append(
  dbGetQuery(
    con,
    'SELECT KHID FROM mesenchyme20hpf LEFT JOIN gene_name ON mesenchyme20hpf.UniqueNAME=gene_name.UniqueNAME'
  ),
  genedat$mesenchyme20hpf
))

mapply(
  dbWriteKey,
  name=names(genedat),
  genedat,
  MoreArgs = list(conn=con, foreign="KH_KY(KHID)")
)

 # bulkRNAseq

# handr <- lrtab('dat/rnaseq/handr/',read.csv,'\\.csv')
# dbWriteGenes(con,'handr_rnaseq',melt.rename(handr,'comparison'))

# ATAC-seq metadata

metadat <- lrtab('data-raw/dat/meta',read.delim,'txt',quote="'")
mapply(dbWriteTable,name=names(metadat),metadat,MoreArgs = list(conn=con,overwrite=T))

ataclib <- dbReadTable(con,'ataclib')

#generate PeakIDs
peakcoord <- paste(
	seqnames(peakome),
	start(peakome)-1,
	end(peakome),sep = ':'
)
peak.id.coord <- setNames(peakome$name,peakcoord)

# ATAC-seq counts
counts <- lrtab(
  'data-raw/dat/counts/',
  read.table
)
dat <- do.call(cbind,sapply(counts,'[',4))
row.names(dat) <- do.call(function(...)paste(...,sep = ':'),counts[[1]][,1:3])

dat <- dat[names(peak.id.coord),]
row.names(dat) <- peak.id.coord#[row.names(dat),]
colnames(dat) <- make.names(sub('_counts.V4','',colnames(dat)))
dat <- as.data.frame(dat)
dat <- dat[,ataclib[,1]]
dbWritePeaks(con,'atacreads',dat,row.names = 'PeakID',overwrite=T)

# add reads per library to metadata
ataclib$reads <- apply(dat,2,sum)
dbWriteTable(con,'ataclib',ataclib,overwrite=T)
