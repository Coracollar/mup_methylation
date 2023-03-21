library(bsseq)
library(DSS)
library(tidyverse)
library(NanoMethViz)
#functions
makeBSseqData <- function (dat, sampleNames) {
  n0 <- length(dat)
  if (missing(sampleNames)) 
    sampleNames <- paste("sample", 1:n0, sep = "")
  alldat <- dat[[1]]
  if (any(alldat[, "N"] < alldat[, "X"], na.rm = TRUE)) 
    stop("Some methylation counts are greater than coverage.\n")
  ix.X <- which(colnames(alldat) == "X")
  ix.N <- which(colnames(alldat) == "N")
  colnames(alldat)[ix.X] <- "X1"
  colnames(alldat)[ix.N] <- "N1"
  if (n0 > 1) {
    for (i in 2:n0) {
      thisdat <- dat[[i]]
      if (any(thisdat[, "N"] < thisdat[, "X"], na.rm = TRUE)) 
        stop("Some methylation counts are greater than coverage.\n")
      ix.X <- which(colnames(thisdat) == "X")
      ix.N <- which(colnames(thisdat) == "N")
      colnames(thisdat)[c(ix.X, ix.N)] <- paste(c("X", 
                                                  "N"), i, sep = "")
      alldat <- merge(alldat, thisdat, all = TRUE, by = c("chr", "pos"))
    }
  }
  alldat <- alldat[order(alldat$chr, alldat$pos), ]
  ix.X <- grep("X", colnames(alldat))
  ix.N <- grep("N", colnames(alldat))
  alldat[is.na(alldat)] <- 0
  M <- as.matrix(alldat[, ix.X, drop = FALSE])
  Cov <- as.matrix(alldat[, ix.N, drop = FALSE])
  colnames(M) <- colnames(Cov) <- sampleNames
  idx <- split(1:length(alldat$chr), alldat$chr)
  M.ordered <- M
  Cov.ordered <- Cov
  pos.ordered <- alldat$pos
  for (i in seq(along = idx)) {
    thisidx = idx[[i]]
    thispos = alldat$pos[thisidx]
    dd = diff(thispos)
    if (min(dd) < 0) {
      warning(paste0("CG positions in chromosome ", names(idx)[i], 
                     " is not ordered. Reorder CG sites.\n"))
      iii = order(thispos)
      M.ordered[thisidx, ] <- M[thisidx, ][iii, ]
      Cov.ordered[thisidx, ] <- Cov[thisidx, ][iii, ]
      pos.ordered[thisidx] <- alldat$pos[thisidx][iii]
    }
  }
  result <- BSseq(chr = alldat$chr, pos = pos.ordered, M = M.ordered, 
                  Cov = Cov.ordered)
  result
}

prepare_for_DSS<-function(query){
  colnames(query) <-c("chr", "pos", "N", "X")
  query$chr <- paste(rep("chr",nrow(query)), query$chr, sep = "")
  query$N<-as.numeric(paste(query$N))
  query$X<-as.numeric(paste(query$X))
  return(query)
}

boxplotting2<-function(BSobj, gr){
  BSobj_gr<-subsetByOverlaps(BSobj, gr)
  methylation_data<-getMeth(BSobj_gr, type="raw")
  methylation_dataframe<-na.omit(as.data.frame(methylation_data))
  methylation_dataframe$Control_mean<-rowMeans(methylation_dataframe[,1:4])  
  methylation_dataframe$CORT_mean<-rowMeans(methylation_dataframe[,5:8])
  methylation_dataframe$Control_sdev<-apply(methylation_dataframe[,1:4],1,sd)  
  methylation_dataframe$CORT_sdev<-apply(methylation_dataframe[,5:8],1,sd) 
  mean(methylation_dataframe$Control_mean)
  mean(methylation_dataframe$CORT_mean)
  methylation_dataframe %>% 
    gather(key="Samples", value="Methylation_frequency") %>%
    ggplot( aes(x=Samples, y=Methylation_frequency, fill=Samples)) +
    geom_boxplot() 
}

#import data
C0<-read.table(file="~/C0/f5c_meth_frec.tsv",  header=T, sep="\t")

C1<-read.table(file="~/C1/f5c_meth_frec.tsv",  header=T, sep="\t")

C2<-read.table(file="~/C2/f5c_meth_frec.tsv", header=T, sep="\t")
  
C3<-read.table(file="~/C3/f5c_meth_frec.tsv",  header=T, sep="\t")

CORT12<-read.table(file="~/CORT12/f5c_meth_frec.tsv",  header=T, sep="\t")

CORT13<-read.table(file="~/CORT13/f5c_meth_frec.tsv", header=T, sep="\t")

CORT14<-read.table(file="~/CORT14/f5c_meth_frec.tsv",  header=T, sep="\t")

CORT15<-read.table(file="~/CORT15/f5c_meth_frec_chr.tsv",  header=T, sep="\t")

C0_DSS<-prepare_for_DSS(C0)
C1_DSS<-prepare_for_DSS(C1)
C2_DSS<-prepare_for_DSS(C2)
C3_DSS<-prepare_for_DSS(C3)
CORT12_DSS<-prepare_for_DSS(CORT12)
CORT13_DSS<-prepare_for_DSS(CORT13)
CORT14_DSS<-prepare_for_DSS(CORT14)
CORT15_DSS<-prepare_for_DSS(CORT15)

rm(C0, C1, C2, C3, CORT12, CORT13, CORT14, CORT15)

BSobjmc <- makeBSseqData( list(C0_DSS, C1_DSS, C2_DSS, C3_DSS, CORT12_DSS, CORT13_DSS, CORT14_DSS, CORT15_DSS),
                        c("C0","C1", "C2", "C3", "CORT12", "CORT13", "CORT14", "CORT15") )

rm(C0_DSS, C1_DSS, C2_DSS, C3_DSS, CORT12_DSS, CORT13_DSS, CORT14_DSS, CORT15_DSS)


#all mup chr4:59,904,830-62,212,385
gr<-GRanges(seqnames = "chr4",
            ranges = IRanges(start = 59904830, end=62212385))


#test significance
dmlTest= DMLtest(BSobjtldr, group1=c("CORT12", "CORT13", "CORT14", "CORT15"), group2=c("C0","C1", "C2", "C3"), smoothing=FALSE)
#find differentially methylated regions
dmrs = callDMR(dmlTest, p.threshold=0.05)

write.table(dmlTest, file= "dmltest_mup.tsv", append = FALSE, sep = "\t", dec = ".", quote = FALSE, row.names = FALSE, col.names = TRUE)
#extract Zhx2 region
gr<-GRanges(seqnames = "chr15",
            ranges = IRanges(start = 57692667, end=57841832))

BSobj_Zhx2<-subsetByOverlaps(BSobj, gr)
dmlTest= DMLtest(BSobj_Zhx2, group1=c("CORT12", "CORT13", "CORT14", "CORT15"), group2=c("C0","C1", "C2", "C3"), smoothing=FALSE)
dmrs = callDMR(dmlTest, p.threshold=0.05) #none found
write.csv(dmlTest, file= "dmltest_Zhx2.csv")

#get methylation frequency matrix
methylation_data<-getMeth(BSobj_mup, type="raw")
rownames(methylation_data)<-paste(BSobj_mup@rowRanges@seqnames, BSobj_mup@rowRanges@ranges@start)

methylation_data2<-getMeth(BSobj_Zhx2, type="raw")
rownames(methylation_data2)<-paste(BSobj_Zhx2@rowRanges@seqnames, BSobj_Zhx2@rowRanges@ranges@start)
alldata<-rbind(methylation_data2, methylation_data)
write.csv(alldata, "mupzhx2_methylation_matrix.csv")

#get methylation frequency matrix
methylation_data<-getMeth(BSobj_mup, type="raw")
rownames(methylation_data)<-paste(BSobj_mup@rowRanges@seqnames, BSobj_mup@rowRanges@ranges@start)
write.table(methylation_data,"~/mup_methylation_matrix.tsv") # table of methylation provided

BSobj_mup<-subsetByOverlaps(BSobjmc, gr)

methylation_data<-getMeth(BSobj_mup, type="raw")
rownames(methylation_data)<-paste(BSobj_mup@rowRanges@seqnames, BSobj_mup@rowRanges@ranges@start)
write.table(methylation_data,"~/mup_methylation_matrixf5c_June22.tsv") # table of methylation provided

methylation_dataframe<-na.omit(as.data.frame(methylation_data))

methylation_dataframe$chr<-as.character(BSobj_mup@rowRanges@seqnames)
methylation_dataframe$start<-BSobj_mup@rowRanges@ranges@start
methylation_dataframe$end<-methylation_dataframe$start+1
methylation_dataframe$Control_mean<-rowMeans(methylation_dataframe[,1:4])  
methylation_dataframe$CORT_mean<-rowMeans(methylation_dataframe[,5:8])

#density plot of Mup region methylation
methylation_dataframe %>% 
  gather(key="Samples", value="Methylation_frequency") %>%
  ggplot( aes(x=Methylation_frequency, colour=Samples)) +
  geom_density()

#boxplot of Mup20 chr4:62,009,410-62,056,143
gr<-GRanges(seqnames = "chr4",
            ranges = IRanges(start = 62009410, end=62056143))
boxplotting2(BSobj_mup, gr)

#Mup20 promoter analysis

#Mup20 promoter 1 chr4	62054106	62054166
#Mup20 promoter 2 chr4	61959968	61960028
gr1<-GRanges("chr4", IRanges(62054106,	62054166))
gr2<-GRanges("chr4", IRanges(61959968,	61960028))
BSobj_p1<-subsetByOverlaps(BSobj_mup, gr1)
p1<-getMeth(BSobj_p1, type="raw")
BSobj_p2<-subsetByOverlaps(BSobj, gr2)
p2<-getMeth(BSobj_p2, type="raw")

df3<-rbind(p1,p2)
df3<-data.frame(df3)
rownames(df3)<-c("p1","p2")
df3<-df3%>%gather(key=Sample, value=Methylation_frequency)
df3$Category<-c("Control","Control","Control","Control","Control","Control","Control","Control",
                "CORT","CORT", "CORT", "CORT","CORT","CORT", "CORT", "CORT")
df3$promoter<-(c("p1", "p2","p1", "p2","p1", "p2","p1", "p2","p1", "p2", "p1", "p2", "p1", "p2", "p1", "p2"))

ggplot(df3, aes(x=promoter, y=Methylation_frequency, fill=Category)) + 
  geom_boxplot(position=position_dodge(0.8))+
  geom_dotplot(binaxis='y', stackdir='center', 
               position=position_dodge(0.8), color="black")

#calculate significance with t.test
m_x<-mean(p1[1:4])
s_x<-sd(p1[1:4])
m_y<-mean(p1[5:8])
s_y<-mean(p1[5:8])
t.test.from.summary.data(m_x, s_x, 4, m_y, s_y, 4)

#plot with NanoMethviZ. Note that the .bgz file was created in advance using the scripts prepare_data_NMV.R and prepare_data_NMV.slurm

methy <-"/data/gpfs/projects/punim1048/RStudio/4vs4/4vs4_methy_data.bgz"

exon_tibble <- get_exons_mus_musculus()

sample <- c("C0", "C1", "C2", "C3","CORT12", "CORT13", "CORT14", "CORT15")
group<-c("Control","Control","Control","Control", "CORT", "CORT", "CORT", "CORT")
sample_anno <- data.frame(sample, group, stringsAsFactors = FALSE)

nmeth_results <- NanoMethResult(methy, sample_anno, exon_tibble)


# plot data

plot_gene(nmeth_results, "Mup2")
plot_gene(nmeth_results, "Mup20")
plot_gene(nmeth_results, "Mup3")
plot_gene(nmeth_results, "Mup18")
plot_gene(nmeth_results, "Mup15")
plot_gene(nmeth_results, "Zhx2")

# whole mup region chr4:59,904,830-62,212,385

plot_region(nmeth_results2, "4", 59904830, 62212385)

