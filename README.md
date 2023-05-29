# mup_methylation

DNA methylation from raw sequencing data was obtained using f5c https://github.com/hasindu2008/f5c. We used a bash script on HPC using 1 GPU and 24 CPUs. 

```bash

# minimap2 --version=2.17-r941
# samtools 1.10 Using htslib 1.10.2
minimap2 -ax map-ont $genomammi basecalls.fastq | samtools view -S -b | samtools sort -o mappings_minimap.sorted.bam &&
samtools index mappings_minimap.sorted.bam
conda deactivate

module load fosscuda
module load f5c/v0.5
export HDF5_PLUGIN_PATH=/home/coracollar/lib/plugins

f5c index -d $fast5folder basecalls.fastq  --iop 24 -t 24

f5c call-methylation --meth-out-version=2 -b $megalodondirectory/mappings_minimap.sorted.bam -g $genomafasta -r basecalls.fastq \
-B 7.0M -K 800 --iop 24 -t 24 > f5c_meth_calls.tsv


f5c meth-freq -i f5c_meth_calls.tsv > f5c_meth_frec.tsv

```

Methylation frequencies obtained from f5c were loaded into R in order to extract the region of interest for differential analysis. We used DSS https://github.com/haowulab/DSS dmltest and findDMR functions to look for differentially methylated regions (no DMR was found in our data for p.threshold 0.05). We used ggplot to plot methylation distributions (boxplot, density plot). Differential methylation at promoters was calculated with a regular ttest.  Functions prepare_for_DSS and boxplotting2 are in the script R_analysis.R.

```R
library(bsseq)
library(DSS)
library(tidyverse)


#mC
C0<-read.table(file="~/C0/20210318_megalodon_4.2.2/f5c_meth_frec.tsv",  header=T, sep="\t")

C1<-read.table(file="~/C1/20210311_megalodon_4.2.2/f5c_meth_frec.tsv",  header=T, sep="\t")

C2<-read.table(file="~/C2_methylation_frec_sinlineachunga.tsv", header=T, sep="\t")
  
C3<-read.table(file="~/C3_f5c_meth_frec.tsv",  header=T, sep="\t")

CORT12<-read.table(file="~/CORT12/megalodon_allfast5_4.2.2/f5c_meth_frec.tsv",  header=T, sep="\t")

CORT13<-read.table(file="~/CORT13/f5c_meth_frec.tsv", header=T, sep="\t")

CORT14<-read.table(file="~/CORT14/20210311_megalodon_4.2.2/f5c_meth_frec.tsv",  header=T, sep="\t")

CORT15<-read.table(file="~/CORT15/20210319_megalodon_4.2.2/f5c_meth_frec_chr.tsv",  header=T, sep="\t")


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

#extract mup region
BSobj_mup<-subsetByOverlaps(BSobjmc, gr)

#test significance
dmlTest= DMLtest(BSobj_mup, group1=c("CORT12", "CORT13", "CORT14", "CORT15"), group2=c("C0","C1", "C2", "C3"), smoothing=FALSE)
#find differentially methylated regions
dmrs = callDMR(dmlTest, p.threshold=0.05) #none found

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

##plot data
#prepare data for ggplot
methylation_dataframe<-na.omit(as.data.frame(methylation_data))
methylation_dataframe$chr<-as.character(BSobj_mup@rowRanges@seqnames)
methylation_dataframe$start<-BSobj_mup@rowRanges@ranges@start
methylation_dataframe$end<-methylation_dataframe$start+1
methylation_dataframe$Control_mean<-rowMeans(methylation_dataframe[,1:4])  
methylation_dataframe$CORT_mean<-rowMeans(methylation_dataframe[,5:8])

#density plot of whole mup region methylation
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

```


Finally, plots of methylation vs genomic position were done using Nanomethviz https://github.com/shians/NanoMethViz, which extracts methylation data from calls derived directly from f5c.

```R
library(NanoMethViz)

#prepare data
C0<-"~/C0/f5c_meth_calls.tsv.gz"

C1<-"~/C1/f5c_meth_calls.tsv.gz"

C2<-"~/C2/f5c_meth_calls.tsv.gz"
  
C3<-"~C3/f5c_meth_calls.tsv.gz"

CORT12<-"~/CORT12/f5c_meth_calls.tsv.gz"

CORT13<-"~/CORT13/f5c_meth_calls.tsv.gz"

CORT14<-"~/CORT14/f5c_meth_calls.tsv.gz"

CORT15<-"~/CORT15/f5c_meth_calls.tsv.gz"


methy_calls <-c(C0, C1, C2, C3, CORT12, CORT13, CORT14, CORT15)

methy_tabix <- file.path("4vs4_methy_data.bgz")
samples <- c("C0", "C1", "C2", "C3","CORT12", "CORT13", "CORT14", "CORT15")
create_tabix_file(methy_calls, methy_tabix, samples)

methy <-"4vs4_methy_data.bgz"

exon_tibble <- get_exons_mus_musculus()
exon_tibble<-exon_tibble[!duplicated(exon_tibble$symbol),]
exon_tibble$chr
exon_tibble$chr<-str_remove(exon_tibble$chr, 'chr')

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

plot_region(nmeth_results, "4", 59904830, 62212385)

#check H19/Igf2 Imprinted region
plot_region(nmeth_results, "7", 142574081,142671474)
gr<-GRanges(seqnames = "chr7",
            ranges = IRanges(start = 142574081, end=142671474))
         
BSobj_imprint<-subsetByOverlaps(BSobjmc, gr)
methylation_data<-getMeth(BSobj_imprint, type="raw")
rownames(methylation_data)<-paste(BSobj_imprint@rowRanges@seqnames, BSobj_imprint@rowRanges@ranges@start)
write.csv(methylation_data,"Igf2_H19_meth_frec.csv" )
dmlTest= DMLtest(BSobj_imprint, group1=c("CORT12", "CORT13", "CORT14", "CORT15"), group2=c("C0","C1", "C2", "C3"), smoothing=FALSE)
write.csv(dmlTest, "dmltest_Igf2H19.csv"
dmrs = callDMR(dmlTest, p.threshold=0.05) #none found
```
