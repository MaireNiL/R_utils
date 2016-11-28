loadRdata <- function(RData, env = new.env()){
  load(RData, env)
  return(env) 
}

dendrogram <- function(){
  indels <- read.table("/Users/ml28/MT_2.0/dendrogram/indels_presence_table.txt",sep="",header=TRUE)
  last_column = (length(indels))
  presence_table = indels[,2:last_column]
  d = dist(t(presence_table)) 
  hc = hclust(d)
  library(rafalib); mypar()
  hcd = as.dendrogram(hc)
  pdf(file='test.pdf',100, 15)
  plot(hcd, cex.lab=0.5)
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

colLab <- function(n) {
  if (is.leaf(n)) {
    a <- attributes(n)
    labCol <- colours[which(labels == a$label)]
    if (is.na(labCol))
      labCol = 8
    attr(n, "nodePar") <- list(a$nodePar, lab.col=labCol, col=labCol, pch=15)
  }
  n
}
  
bafPlot <- function() {
for (chrom in chromList) {
    chrom.idx = snvs.metadata$CHROM == chrom
    for (i in 1:ncol(snvs.nr)) {
        png(paste0("Chrom", chrom, "_", samples[i], "_BAFplots_AllSamples_tonly_test.png"), width=6000,height=1500)
        chrom.tonly.idx = tumour.only.idx[chrom.idx]
        chrom.snponly.idx = snps.idx[chrom.idx]
        chrom.tunique.idx = tumour.unique.idx[chrom.idx]
        vaf=snvs.vaf[,i]
        plot((snvs.metadata$POS)[chrom.idx][chrom.snponly.idx], vaf[chrom.idx][chrom.snponly.idx],main=paste0("Chromosome ", chrom, ", ", samples[i]), xlab="SNV index", ylab="BAF", pch=20, col="black", cex=1.7, cex.main=1.6, cex.axis=1.4, cex.lab=1.4)
        points(snvs.metadata$POS[chrom.idx][chrom.tonly.idx], vaf[chrom.idx][chrom.tonly.idx], pch=20, col="red", cex=1.7)
        points(snvs.metadata$POS[chrom.idx][chrom.tunique.idx], vaf[chrom.idx][chrom.tunique.idx], pch=20, col="chartreuse", cex=2.5)
        par(xpd=TRUE)
        legend(y=1.5, x=0, horiz=T,  legend=c("SNPS","Tumour-only", "Tumour-unique"), pch=c(20,20,20), col=c("black","red","green"), pt.cex=2, cex=1.6, bty="n", xpd=T)
        dev.off()
    }
}
 
dStat <- function(dstatsInput, pdfOutput) {
library(plotly)
df <- read.table(dstatsInput, header=FALSE)
pdf(file=pdfOutput, 10, 100, useDingbats = FALSE)
ggplot(df, aes(V5, V6)) +
  stat_boxplot(geom ='errorbar') + 
  geom_boxplot(position=position_dodge(width=0.8)) +
  coord_flip() +
  ggtitle("D-statistics D(Outgroup, Admix Test Pop, Pop X, Pop Y)") +
  xlab("Pop X") + 
  ylab("D") +
  geom_hline(yintercept=0.0, colour="red", linetype="dashed", size=1 )

# Note that plotly will automatically add horozontal lines to the whiskers
# if we need to add colours use this
ggplot(df, aes(V5, V6), fill=V6) +
  geom_boxplot(position=position_dodge(width=0/8)) +
  scale_color_gradientn(colours = rainbow(276)) +
  coord_flip() +
  geom_hline(yintercept=0.0, colour="red", linetype="dashed", size=1 )

dev.off()
}

logRplot <-function() {
for (group in chromList) {
    chrom.idx = snvs.metadata$CHROM == group
    sample_pairs <- cbind(samples, "xxxH-Dog") # where xxH-Dog is the host we want to normalize to
    for (i in 1:nrow(sample_pairs)) {
        par(mfrow=c(1,2)) 
        # update this using layout
        sample = cbind(snvs.nr[,sample_pairs[i,1]], snvs.nr[,sample_pairs[i,2]])
        logR = apply(sample, 1, function(dat) {
        nT = dat[1]
        nH = dat[2]
        log2(nT / nH)
        })
     plot(snvs.metadata$POS[chrom.idx], logR[chrom.idx],yaxt="n", pch=20, col="gray45", main=paste0("BAF and logR - ", sample_pairs[i,1]), cex=1.3, cex.main=2, cex.axis=2,xlab="Genomic coordinate", ylab="logR")
      abline(h=median(na.omit(logR)))     
    }
    dev.off()
}}

function <- bafplotOriginal() {
# Code from AOB 18/08/2016
# BAFPLOTS PER CHROMOSOME USING CURRENT DATA OBJECTS 
# Load data
# ......snvs.metadata, snvs.nr, snvs.nv, snvs.vaf

# Create output directory
outdir = paste0("~/Desktop/", Sys.Date(), "_BAFplots")
dir.create(outdir, showWarnings=F)

# A) ONE FILE PER SAMPLE, ONE PAGE PER CHROMOSOME
##################################################

for (i in 1:ncol(snvs.nr)) {
    pdf(paste0(outdir, "/", samples[i], "_BAFplots_AllChrom.pdf"), width=40, height=10)
    par(mar=c(5.1, 4.9, 5, 2.1))
    
    # Take all sample VAFs
    vaf = snvs.vaf[,i]
    
    for (chrom in c(1:38, "X")) {
        # Take variants in the chromosome
        chrom.idx = snvs.metadata$CHROM == chrom
        
        # First, plot all variants
        plot(1:sum(chrom.idx), vaf[chrom.idx], main=paste0(samples[i], ", chomosome ", chrom), xlab="SNV index", ylab="BAF", pch=20, 
             col="gray72", cex=1.7, cex.main=1.6, cex.axis=1.4, cex.lab=1.4)
        
        # Second, plot tumour-only variants in black
        chrom.tonly.idx = tumour.only.idx[chrom.idx]
        points((1:sum(chrom.idx))[chrom.tonly.idx], vaf[chrom.idx][chrom.tonly.idx], pch=18, col="black", cex=2.2)
        
        # Add legend
        legend(y=1.15, x=0, horiz=T,  legend=c("Tumour-only     ", "Host+tumour"), pch=c(18,20), col=c("black", "grey"), pt.cex=2, cex=1.6, bty="n", xpd=T)
    }
    dev.off()
}

# B) ONE FILE PER CHROMOSOME, ONE PAGE PER SAMPLE
##################################################

for (chrom in 1) { #c(1:38, "X")) {
    pdf(paste0(outdir, "/Chrom", chrom, "_BAFplots_AllSamples.pdf"), width=40, height=10)
    par(mar=c(5.1, 4.9, 5, 2.1))
    
    # Take variants in the chromosome
    chrom.idx = snvs.metadata$CHROM == chrom
    
    for (i in 1:ncol(snvs.nr)) {
        # Take all sample VAFs
        vaf = snvs.vaf[,i]
        
        # First, plot all variants
        plot(1:sum(chrom.idx), vaf[chrom.idx], main=paste0("Chromosome ", chrom, ", ", samples[i]), xlab="SNV index", ylab="BAF", pch=20, 
             col="gray72", cex=1.7, cex.main=1.6, cex.axis=1.4, cex.lab=1.4)
        
        # Second, plot tumour-only variants in black
        chrom.tonly.idx = tumour.only.idx[chrom.idx]
        points((1:sum(chrom.idx))[chrom.tonly.idx], vaf[chrom.idx][chrom.tonly.idx], pch=18, col="black", cex=2.2)
        
        # Add legend
        legend(y=1.15, x=0, horiz=T,  legend=c("Tumour-only     ", "Host+tumour"), pch=c(18,20), col=c("black", "grey"), pt.cex=2, cex=1.6, bty="n", xpd=T)
    }
    
    dev.off()
}}
 
