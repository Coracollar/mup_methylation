library(NanoMethViz)

#prepare data
C0<-"~/C0_160921_megalodon_5mc_5hmc/per_read_modified_base_calls.db.gz"

C1<-"~/C1_1011921_megalodon_5mc_5hmc/per_read_modified_base_calls.db.gz"

C2<-"~/C2_181121_megalodon_5mc_5hmc/per_read_modified_base_calls.db.gz"
  
C3<-"~C3_261121_megalodon_5mc_5hmc/per_read_modified_base_calls.db.gz"

CORT12<-"~/CORT12_241121_megalodon_5mc_5hmc/per_read_modified_base_calls.db.gz"

CORT13<-"~/CORT13_281021_megalodon_5mc_5hmc/per_read_modified_base_calls.db.gz"

CORT14<-"~/CORT14_251121_megalodon_5mc_5hmc/per_read_modified_base_calls.db.gz"

CORT15<-"~/CORT15_261021_megalodon_5mc_5hmc/per_read_modified_base_calls.db.gz"


methy_calls <-c(C0, C1, C2, C3, CORT12, CORT13, CORT14, CORT15)

methy_tabix <- file.path("megalodon4vs4_methy_data.bgz")
samples <- c("C0", "C1", "C2", "C3","CORT12", "CORT13", "CORT14", "CORT15")
create_tabix_file(methy_calls, methy_tabix, samples)

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
