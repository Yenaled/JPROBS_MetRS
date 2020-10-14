# analyze.r
# Input file should be supplied as data.xlsx in the current working directory

require(openxlsx)
require(limma)
require(imputeLCMD)

# Set thresholds
lfcThresh <<- 1
fdr <<- 0.05

# Make output folder if it doesn't already exist

output_dir <<- "./output/"
dir.create(output_dir, showWarnings = FALSE)

# Which column contains the Scaffold DIA program's computed log2 fold changes
lfc_column <<- "log2(FC)2"

# Function to impute missing values
# Adopted from https://www.rdocumentation.org/packages/imputeLCMD/versions/2.0/topics/impute.QRILC

set.seed(42)

imputeMissingValues <- function(data, nFeatures) {
  tune.sigma = 1
  QR.obj = list()
  curr.sample = data
  pNAs = length(which(is.na(curr.sample)))/length(curr.sample)
  upper.q = 0.95
  q.normal = qnorm(seq(pNAs, upper.q, (upper.q - pNAs)/(upper.q * 10000)), mean = 0, sd = 1)
  q.curr.sample = quantile(curr.sample, probs = seq(0, upper.q, 1e-04), na.rm = T)
  temp.QR = lm(q.curr.sample ~ q.normal)
  mean.CDD = temp.QR$coefficients[1]
  sd.CDD = as.numeric(temp.QR$coefficients[2])
  data.to.imp = rtmvnorm(n = nFeatures, mean = mean.CDD, 
                           sigma = sd.CDD * tune.sigma, upper =
                           qnorm(pNAs, mean = mean.CDD, sd = sd.CDD), algorithm = c("gibbs"))
    curr.sample.imputed = curr.sample
    curr.sample.imputed[which(is.na(curr.sample))] = data.to.imp[which(is.na(curr.sample))]
    return(curr.sample.imputed)
}

# Function to clean the data loaded from the spreadsheet

cleanSpreadsheetData <- function(data, group1=NULL, group2=NULL) {
  data <- data[!is.na(data[,"Accession Number"]),]
  # Replace NQ with NC in the log2 fold change column
  data[data[,lfc_column] == "NQ",lfc_column] <- "NC"
  # Replace "Reference Missing" with NC in the log2 fold change column
  data[data[,lfc_column] == "Reference Missing",lfc_column] <- "NC"
  if (!is.null(group1) && !is.null(group2)) {
    # If >=50% of values are missing in a group, fix the log2 fold change column if it's negative (because parantheses got converted into minus sign)
    group1_fix <- rowSums(data[,group1] == "NQ") >= length(group1)/2 & rowSums(data[,group1] == "NQ") < length(group1)-1
    group1_fix <- names(group1_fix[group1_fix == TRUE])
    group1_fix <- rownames(data) %in% group1_fix & !grepl("\\(", data[,lfc_column]) & data[,lfc_column] != "NC"
    data[group1_fix,lfc_column] <- abs(as.numeric(data[group1_fix,lfc_column]))
    group2_fix <- rowSums(data[,group2] == "NQ") >= length(group2)/2 & rowSums(data[,group2] == "NQ") < length(group2)-1
    group2_fix <- names(group2_fix[group2_fix == TRUE])
    group2_fix <- rownames(data) %in% group2_fix & !grepl("\\(", data[,lfc_column]) & data[,lfc_column] != "NC"
    data[group2_fix,lfc_column] <- abs(as.numeric(data[group2_fix,lfc_column]))
  }
  # Remove parantheses for numbers that are enclosed within parantheses in the log2 fold change column
  data[grepl("\\(", data[,lfc_column]),lfc_column] <- gsub("\\)", "", gsub("\\(", "", data[grepl("\\(", data[,lfc_column]),lfc_column]))
  data[data[,lfc_column] != "NC",lfc_column] <- sprintf("%.2f",round(as.numeric(data[data[,lfc_column] != "NC",lfc_column]), 2))
  return(data)
}

# Load boncat enriched data

e_data <- read.xlsx("data.xlsx", check.names=FALSE, sep.names=" ", startRow=13, sheet=1)
e_met <- c("31CE","30CE","32CE","33CE")
e_anl <- c("40AE","ANL1E","38AE","29AE")
e_data <- cleanSpreadsheetData(e_data, e_met, e_anl)

# Load whole bulk 

w_data <- read.xlsx("data.xlsx", check.names=FALSE, sep.names=" ", startRow=13, sheet=2)
w_met <- c("31CW", "30CW", "32CW", "33CW")
w_anl <- c("40AW", "ANL1W", "38AW", "29AW")
w_data <- cleanSpreadsheetData(w_data, w_met, w_anl)

# Note: e_met and w_met are paired samples and are ordered to reflect the pairings; same thing with e_anl and w_anl

# Prepare data for limma analysis

e_data_formatted <- e_data[,c("Accession Number", e_met, e_anl)]
rownames(e_data_formatted) <- e_data_formatted[,1]
e_data_formatted[,1] <- NULL
e_data_formatted[] <- data.frame(sapply(e_data_formatted, function(x) as.numeric(x)))

w_data_formatted <- w_data[,c("Accession Number", w_met, w_anl)]
rownames(w_data_formatted) <- w_data_formatted[,1]
w_data_formatted[,1] <- NULL
w_data_formatted[] <- data.frame(sapply(w_data_formatted, function(x) as.numeric(x)))

data_combined <- merge(e_data_formatted, w_data_formatted, by="row.names", all=FALSE)
rownames(data_combined) <- data_combined[,1]
data_combined[,1] <- NULL
data_combined <- data_combined[,c(w_anl, e_anl)]

# Transform (log10) data

e_data_normalized <- log10(e_data_formatted)
w_data_normalized <- log10(w_data_formatted)
data_normalized <- log10(data_combined)

# Volcano plot function

volcanoplot <- function (res, lfcthresh=1, sigthresh=0.05, highlight=NULL, main="Volcano Plot", legendpos="bottomright", textcx=1, ...) {
  res <- res[,c("logFC", "P.Value", "adj.P.Val", "gene_symbol")]
  colnames(res) <- c("log2FoldChange", "pvalue", "padj", "Gene")
  #res$pvalue <- res$padj
  ymax <- 10
  res[!is.na(res$pvalue) & res$pvalue < 10^-ymax, "pvalue"] <- 10^-ymax
  scale_axis <- c(-12, -8, -4, 0, 4, 8, 12)
  xlim <- c(-12, 12)
  if (max(abs(res$log2FoldChange), na.rm=TRUE) < 6) {
    xlim <- c(-6, 6)
    scale_axis <- c(-6, -3, 0, 3, 6)
  }
  with(res, plot(log2FoldChange, -log10(pvalue), xaxt='n', pch=20, cex.lab=1.3, cex.axis=1.7, xlab="Log2 Fold Change", ylab="-Log10(p-value)", col="black", xlim=xlim, ylim=c(0,ymax), main=main, ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="blue", ...))
  if (!is.null(highlight)) {
    with(subset(res, Gene %in% highlight), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
    print(main)
    highlight_res <- res[res$Gene %in% highlight,]
    highlight_res <- highlight_res[order(highlight_res$log2FoldChange),]
    rownames(highlight_res) <- 1:nrow(highlight_res)
    print(highlight_res)
  }
  #with(res, abline(lty=2, h=-log10(sigthresh)))
  highest_sig_p <- max(res[!is.na(res$padj) & res$padj < sigthresh,"pvalue"])
  with(res, abline(lty=2, h=-log10(highest_sig_p)))
  with(res, abline(lty=2, v=-lfcthresh))
  with(res, abline(lty=2, v=lfcthresh))
  with(res, axis(side=1,cex.lab=1.5, cex.axis=2, at=scale_axis))
  print(paste(main, " UP: ", nrow(res[!is.na(res$padj) & !is.na(res$log2FoldChange) & res$padj < sigthresh & res$log2FoldChange > lfcthresh,]), sep=""))
  print(paste(main, " DOWN: ", nrow(res[!is.na(res$padj) & !is.na(res$log2FoldChange) & res$padj < sigthresh & res$log2FoldChange < -lfcthresh,]), sep=""))
}

# Function to merge results with gene symbols
# Note: Mapping from protein accessions to gene symbols was done via https://www.uniprot.org/uploadlists/

accession_to_gene_mapping <<- read.table("protein_accessions_to_genes.tab.txt", stringsAsFactors=FALSE, sep="\t", header=TRUE)
colnames(accession_to_gene_mapping) <- c("protein", "gene_symbol")
accession_to_gene_mapping <- accession_to_gene_mapping[!duplicated(accession_to_gene_mapping$protein),] # Remove one-to-many protein-gene mappings (e.g. histones)
combine_results <- function(results, original_data) {
  results$full_accession <- rownames(results)
  results$accession <- gsub(".*\\|(.+)\\|.*", "\\1", rownames(results))
  results <- merge(results, accession_to_gene_mapping, by.x="accession", by.y="protein", all.x=TRUE, all.y=FALSE)
  rownames(results) <- results[,"full_accession"]
  results[,"full_accession"] <- NULL
  results$accession <- NULL
  if (!is.null(original_data)) {
    original_order <- original_data[,"Accession Number"]
    intersect_cols <- intersect(colnames(results), colnames(original_data))
    original_data[,intersect_cols] <- NULL
    results <- merge(original_data, results, all=TRUE, by.x="Accession Number", by.y="row.names")
    results <- results[match(original_order, results[,"Accession Number"]),] # Re-order rows
    results <- results[,unique(c(colnames(original_data), colnames(results)))] # Re-order columns
  } else {
    results[,"Accession Number"] <- rownames(results)
    results <- results[,colnames(results)[c(ncol(results), 1:(ncol(results)-1))]]
    results <- results[order(results[,"adj.P.Val"]),]
  }
  return(results)
}

# Function to do differential analysis

doDiff <- function(data, original_data, title, conditions, wb, highlight=NULL, paired=FALSE) {
  fname <- gsub(" ", "_", title)
  columns <- unique(c(paste("t", conditions, sep="")))
  if (paired) {
    num_samples <- length(conditions) / 2
    pairings <- factor(c(1:num_samples, 1:num_samples))
    design <- model.matrix(~0+conditions+pairings)
    colnames(design)[1:length(columns)] <- columns
  } else {
    design <- model.matrix(~0+conditions)
    colnames(design) <- columns
  }
  contrast.list <- paste(columns[2], columns[1], sep="-")
  contrast <- do.call(makeContrasts, c(contrast.list, list(levels=design)))
  fit <- lmFit(as.matrix(data), design)
  fit2 <- contrasts.fit(fit, contrast)
  fit2 <- eBayes(fit2)
  results <- topTable(fit2,adjust.method="BH",n=100000,p.value=1, coef=1)
  results <- combine_results(results, original_data)
  condition2 <- as.numeric(rownames(design[design[,columns[2]] == 1,]))
  condition1 <- as.numeric(rownames(design[design[,columns[1]] == 1,]))
  
  # Compute median log2 fold changes
  data_median1 <- apply(data[,condition1], 1, FUN = median, na.rm=TRUE)
  data_median2 <- apply(data[,condition2], 1, FUN = median, na.rm=TRUE)
  data$logFC <- (data_median2 - data_median1) / log10(2)
  ### If all samples or all-but-one-sample in a group is NA, then make logFC NA
  group1_fix <- rowSums(is.na(data[,condition1])) >= length(condition1)-1
  group1_fix <- names(group1_fix[group1_fix == TRUE])
  data[group1_fix,"logFC"] <- NA
  group2_fix <- rowSums(is.na(data[,condition2])) >= length(condition2)-1
  group2_fix <- names(group2_fix[group2_fix == TRUE])
  data[group2_fix,"logFC"] <- NA
  ### If >=50% of values are missing in a group, impute missing values for the log2 fold changes
  group1_fix2 <- rowSums(is.na(data[,condition1])) >= length(condition1)/2 & rowSums(is.na(data[,condition1])) < length(condition1)-1
  group1_fix2 <- names(group1_fix2[group1_fix2 == TRUE])
  group2_fix2 <- rowSums(is.na(data[,condition2])) >= length(condition2)/2 & rowSums(is.na(data[,condition2])) < length(condition2)-1
  group2_fix2 <- names(group2_fix2[group2_fix2 == TRUE])
  group1_fix2 <- group1_fix2[!(group1_fix2 %in% c(group1_fix, group2_fix))]
  group2_fix2 <- group2_fix2[!(group2_fix2 %in% c(group1_fix, group2_fix))]
  data_to_impute <- data[unique(c(group1_fix2,group2_fix2)),]
  nFeatures = length(condition1) + length(condition2)
  for (i in rownames(data_to_impute)) {
    if (i %in% group1_fix2) {
      data_to_impute[i,condition1] <- imputeMissingValues(data_to_impute[i,condition1], nFeatures)
    }
    if (i %in% group2_fix2) {
      data_to_impute[i,condition2] <- imputeMissingValues(data_to_impute[i,condition2], nFeatures)
    }
    xx1 <- as.numeric(data_to_impute[i,condition1])
    xx2 <- as.numeric(data_to_impute[i,condition2])
    xx1[is.na(xx1)] <- min(0, xx1[!is.na(xx1)])
    xx2[is.na(xx2)] <- min(0, xx2[!is.na(xx2)])
    data_median1 <- median(xx1)
    data_median2 <- median(xx2)
    data[i,"logFC"] <- ((data_median2 - data_median1) / log10(2))
  }
  ### Fix instances where there's <50% of values missing in a group
  group1_fix3 <- rowSums(is.na(data[,condition1])) < length(condition1)/2 & rowSums(is.na(data[,condition1])) > 0
  group1_fix3 <- names(group1_fix3[group1_fix3 == TRUE])
  group2_fix3 <- rowSums(is.na(data[,condition2])) < length(condition2)/2 & rowSums(is.na(data[,condition2])) > 0
  group2_fix3 <- names(group2_fix3[group2_fix3 == TRUE])
  group1_fix3 <- group1_fix3[!(group1_fix3 %in% c(group1_fix, group2_fix, group1_fix2, group2_fix2))]
  group2_fix3 <- group2_fix3[!(group2_fix3 %in% c(group1_fix, group2_fix, group1_fix2, group2_fix2))]
  data_to_impute <- data[unique(c(group1_fix3,group2_fix3)),]
  for (i in rownames(data_to_impute)) {
    xx1 <- as.numeric(data_to_impute[i,condition1])
    xx2 <- as.numeric(data_to_impute[i,condition2])
    xx1[is.na(xx1)] <- min(0, xx1[!is.na(xx1)])
    xx2[is.na(xx2)] <- min(0, xx2[!is.na(xx2)])
    data_median1 <- median(xx1)
    data_median2 <- median(xx2)
    data[i,"logFC"] <- ((data_median2 - data_median1) / log10(2))
  }
  # End computing fold changes
  
  data[,"Accession Number"] <- rownames(data)
  colnames(results)[colnames(results) == 'logFC'] <- 'beta'
  results <- merge(results, data[,c("Accession Number", "logFC")], by="Accession Number")
  #print(results[1:5,])
  #print(sum(is.na(results[,"logFC"])))
  #print(nrow(results[results[,lfc_column] == "NC",]))
  #pdf("test.pdf")
  #print((results[!is.na(as.numeric(results[,lfc_column])) & as.numeric(results[,lfc_column]) > 25,]))
  #print((results[!is.na(as.numeric(results[,lfc_column])) & as.numeric(results[,lfc_column]) < -0.25 & as.numeric(results[,"logFC"]) > 0.1,]))
  #plot(as.numeric(results[,lfc_column]), as.numeric(results$logFC))
  #print(cor.test(as.numeric(results[,lfc_column]), as.numeric(results$logFC)))
  #dev.off()
  

  # Sort by original data
  if (!is.null(original_data)) {
    original_order <- original_data[,"Accession Number"]
    results <- results[match(original_order, results[,"Accession Number"]),] # Re-order rows
    results[,c(1,2)] <- results[,c(2,1)]
  }

  addWorksheet(wb, title)
  writeData(wb, sheet = title, results, rowNames = FALSE)
  if (lfc_column %in% colnames(results)) {
    results$logFC <- as.numeric(results[,lfc_column])
  }
  pdf(paste(output_dir, "volcano_", fname, ".pdf", sep=""))
  volcanoplot(results, lfcthresh=lfcThresh, sigthresh=fdr, main=title, highlight=highlight)
  device <- dev.off()
  up_genes <- results[!is.na(results[,"adj.P.Val"]) & results[,"adj.P.Val"] < fdr & !is.na(results$logFC) & results$logFC > lfcThresh,"gene_symbol"]
  down_genes <- results[!is.na(results[,"adj.P.Val"]) & results[,"adj.P.Val"] < fdr & !is.na(results$logFC) & results$logFC < -lfcThresh,"gene_symbol"]
  writeLines(up_genes, paste(output_dir, "UPgenes_", fname, ".txt", sep=""))
  writeLines(down_genes, paste(output_dir, "DOWNgenes_", fname, ".txt", sep=""))
}
wb <- createWorkbook("Supplementary Table")
sink(paste(output_dir, "highlights.txt", sep=""), split=TRUE)

# Analyze BONCAT ANL vs. BONCAT MET

data <- e_data_normalized
original_data <- e_data
title <- "BONCAT-ANL vs BONCAT-Met"
conditions <- factor(c(rep(0, length(e_met)), rep(1, length(e_anl))))
highlight <- c("Kras", "Hmgb2", "Yap1", "Lgals3", "Hmgb1") # Note: Lgals3 = LEG3
doDiff(data, original_data, title, conditions, wb, highlight)

# Analyze BULK ANL vs. BULK MET

data <- w_data_normalized
original_data <- w_data
title <- "Bulk-ANL vs Bulk-Met"
conditions <- factor(c(rep(0, length(w_met)), rep(1, length(w_anl))))
doDiff(data, original_data, title, conditions, wb, highlight)

# Analyze BONCAT ANL vs. BULK ANL

data <- data_normalized
original_data <- NULL
title <- "BONCAT-ANL vs Bulk-ANL"
conditions <- factor(c(rep(0, length(w_anl)), rep(1, length(e_anl))))
highlight <- c("Kras", "Mki67") # Note: Mki67 = KI67
doDiff(data, original_data, title, conditions, wb, highlight, paired=TRUE)

# Save analysis to Excel

sink()
saveWorkbook(wb, paste(output_dir, "analysis.xlsx", sep=""), overwrite = TRUE)

# Print some information about the boncat ANL vs. MET enrichment

results <- read.xlsx(paste(output_dir, "analysis.xlsx", sep=""), check.names=FALSE, sep.names=" ", sheet=1)
lfc_column <- "logFC"
results[,"adj.P.Val"] <- as.numeric(results[,"adj.P.Val"])
results[,lfc_column] <- as.numeric(results[,lfc_column])
up <- results[!is.na(results[,"adj.P.Val"]) & !is.na(results[,lfc_column]) & results[,"adj.P.Val"] < fdr & results[,lfc_column] > lfcThresh,"Protein Name"]
down <- results[!is.na(results[,"adj.P.Val"]) & !is.na(results[,lfc_column]) & results[,"adj.P.Val"] < fdr & results[,lfc_column] < -lfcThresh,"Protein Name"]

missing_value_label <<- "NQ"
aa <- (e_data[apply(e_data[,c(e_met)], 1, function(x) { !(length(unique(x)) == 1 & x[1] == missing_value_label) } ),"Protein Name"])
bb <- (e_data[apply(e_data[,c(e_anl)], 1, function(x) { !(length(unique(x)) == 1 & x[1] == missing_value_label) } ),"Protein Name"])
cc <- (e_data[apply(e_data[,c(e_met,e_anl)], 1, function(x) { !(length(unique(x)) == 1 & x[1] == missing_value_label) } ),"Protein Name"])

all_proteins <- unique(c(aa,bb))
proteins_met_only <- unique(c(down, aa[!(aa %in% intersect(aa,bb))]))
proteins_anl_only <- unique(c(up, bb[!(bb %in% intersect(aa,bb))]))
#proteins_met_only <- unique(c(aa[!(aa %in% intersect(aa,bb))]))
#proteins_anl_only <- unique(c(bb[!(bb %in% intersect(aa,bb))]))
if (!identical(cc[order(cc)], all_proteins[order(all_proteins)])) {
    stop("Error1")
}
if (length(intersect(proteins_met_only, proteins_anl_only)) != 0) {
    stop("Error2")
}
proteins_met_anl <- unique(c(proteins_met_only, proteins_anl_only))
proteins_intersection <- all_proteins[!(all_proteins %in% proteins_met_anl)]
intersection_size <- length(all_proteins) - length(proteins_met_only) - length(proteins_anl_only)
if (length(proteins_intersection) != intersection_size) {
    stop("Error3")
}
print(paste("E-MET PROTEINS ENRICHED:", length(proteins_met_only))) # 35
print(paste("E-ANL PROTEINS ENRICHED:", length(proteins_anl_only))) # 3727
print(paste("INTERSECTION:", intersection_size)) # 645
print(paste("Total number of proteins:", length(all_proteins))) # 4407
print("----")
print(paste("E-MET PROTEINS IDENTIFIED:", length(aa))) # 4386
print(paste("E-ANL PROTEINS IDENTIFIED:", length(bb))) # 4405

