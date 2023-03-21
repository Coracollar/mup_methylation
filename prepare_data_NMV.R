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


