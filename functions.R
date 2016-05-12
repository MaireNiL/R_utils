dendrogram <- function(){
  indels<-read.table("/Users/ml28/MT_2.0/dendrogram/indels_presence_table.txt",sep="",header=TRUE)
  last_column=(length(indels))
  presence_table=indels[,2:last_column]
  d = dist(t(presence_table)) 
  hc = hclust(d)
  library(rafalib); mypar()
  hcd = as.dendrogram(hc)
  pdf(file='test.pdf',100, 15)
  plot(hcd,cex.lab=0.5)
dev.off()
}

strand.bias <- function(VCF_Input, variant_support_edge, min_reads, min_percentage){
  # Adapted from code by Maximilian Stammnitz, mrs72@cam.ac.uk
  # 1. Build 3 columns: N/TC, NF/N, NR/N from Info-field
  
  cat("\n Isolating strand-Coverages...")
  N   <- str_split_fixed(as.character(VCF_Input[,"INFO"]), ";", 20)[,18]
  NF  <- str_split_fixed(as.character(VCF_Input[,"INFO"]), ";", 20)[,8]
  NR  <- str_split_fixed(as.character(VCF_Input[,"INFO"]), ";", 20)[,9]
  TC  <- str_split_fixed(as.character(VCF_Input[,"INFO"]), ";", 20)[,15]
  N_TC <- paste(str_split_fixed(N, "=", 2)[,2], str_split_fixed(TC, "=", 2)[,2], sep="/")
  NF_N <- paste(str_split_fixed(NF, "=", 2)[,2], str_split_fixed(N, "=", 2)[,2], sep="/")
  NR_N <- paste(str_split_fixed(NR, "=", 2)[,2], str_split_fixed(N, "=", 2)[,2], sep="/")
  names <- colnames(VCF_Input)
  VCF_Input <- cbind(VCF_Input, N_TC, NF_N, NR_N)
  colnames(VCF_Input) <- c(names, "Support Coverage", "Forward Support", "Reverse Support")
  
  # 2. Build two position vectors: which total support is ≤ or > variant_support_edge
  cat("\n Filtering...")
  above <- which(suppressWarnings(as.numeric(str_split_fixed(VCF_Input[,"Support Coverage"], "/", 2)[,1]))>variant_support_edge)
  below_or_equal <- which(suppressWarnings(as.numeric(str_split_fixed(VCF_Input[,"Support Coverage"], "/", 2)[,1]))<=variant_support_edge)
  
  # 3. Apply thresholds to positions in either vector, but on the whole set!
  
  # above: which Forward Support or Reverse Support are < 20% ?
  above.fw <- as.numeric(str_split_fixed(as.character(VCF_Input[above,"Forward Support"]),"/",2)[,1])/as.numeric(str_split_fixed(as.character(VCF_Input[above,"Forward Support"]),"/",2)[,2])
  above.rv <- as.numeric(str_split_fixed(as.character(VCF_Input[above,"Reverse Support"]),"/",2)[,1])/as.numeric(str_split_fixed(as.character(VCF_Input[above,"Reverse Support"]),"/",2)[,2])
  above <- above[which(above.fw<min_percentage | above.rv<min_percentage)]
  
  ## below: which Forward Support or Reverse Support are < 2 reads ?
  below_or_equal.fw <- as.numeric(str_split_fixed(as.character(VCF_Input[below_or_equal,"Forward Support"]),"/",2)[,1])
  below_or_equal.rv <- as.numeric(str_split_fixed(as.character(VCF_Input[below_or_equal,"Reverse Support"]),"/",2)[,1])
  below_or_equal <- below_or_equal[which(below_or_equal.fw<min_reads | below_or_equal.rv<min_reads)]
 
  # variants which pass the tresholds
  fail <- sort(union(above,below_or_equal))
  fail_10 <- sort(below_or_equal)  
  # 4. Count the numbers of removed variants (and percentage of the total calls)
  cat("\n \n Variants removed due to percentage underpassing: ", length(which(above.fw<min_percentage | above.rv<min_percentage)), "of", length(N))
  cat("\n Variants removed due to read underpassing: ", length(which(below_or_equal.fw<min_reads | below_or_equal.rv<min_reads)), "of", length(N))
  cat("\n Total removed: ", round(c(length(which(above.fw<min_percentage | above.rv<min_percentage))+length(which(below_or_equal.fw<min_reads | below_or_equal.rv<min_reads)))/length(N),4)*100, "%")
  
  # 5. Output
 VCF_Input <- VCF_Input[-fail,]
 return(VCF_Input)
}
