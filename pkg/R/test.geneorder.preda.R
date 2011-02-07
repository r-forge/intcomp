test.geneorder.preda <- function(ge, cn, Labels, nperm, ge.qval.threshold=0.05, cn.qval.threshold=0.01, smoothMethod="spline", ge.smoothStatistic.threshold.up=0.5, ge.smoothStatistic.threshold.down=-0.5, cn.smoothStatistic.threshold.gain=0.1, cn.smoothStatistic.threshold.loss=-0.1, chromosomes=1:22, cancerGenes){
require(PREDA)

# prepare annotation for gene expression
GEdataFile <- as.data.frame(cbind(rownames(ge$data),ge$data))
GEStatisticsForPREDA <- StatisticsForPREDAFromdataframe(StatisticsForPREDA_dataframe=GEdataFile,ids_column=1,statistic_columns=2:(ncol(ge$data)+1), testedTail="both")

GEGenomicsAnnotations <- new("GenomicAnnotations", 
        ids = as.character(rownames(ge$data)), 
        chr = as.numeric(ge$info$chr), 
        start = as.numeric(ge$info$start), 
        end = as.numeric(ge$info$end),
        strand = as.numeric(rep(1,nrow(ge$data))), 
        chromosomesNumbers = chromosomes,
        chromosomesLabels = as.character(chromosomes))
        
GEGenomicsAnnotationsForPREDA <- GenomicAnnotations2GenomicAnnotationsForPREDA(GEGenomicsAnnotations, reference_position_type="median")

SODEGIRGEDataForPREDA <- MergeStatisticAnnotations2DataForPREDA(GEStatisticsForPREDA, GEGenomicsAnnotationsForPREDA, sortAndCleanNA = TRUE)

# gene expression analysis
SODEGIRGEanalysisResults <- PREDA_main(SODEGIRGEDataForPREDA, smoothMethod = smoothMethod, nperms = nperm)

SODEGIR_GE_UP <- PREDAResults2GenomicRegions(SODEGIRGEanalysisResults,
   qval.threshold = ge.qval.threshold, smoothStatistic.tail = "upper",
   smoothStatistic.threshold = ge.smoothStatistic.threshold.up)
SODEGIR_GE_DOWN <- PREDAResults2GenomicRegions(SODEGIRGEanalysisResults,
   qval.threshold = ge.qval.threshold, smoothStatistic.tail = "lower",
   smoothStatistic.threshold = ge.smoothStatistic.threshold.down)

# prepare annotation for copy number
CNdataFile <- as.data.frame(cbind(rownames(cn$data),cn$data))
CNStatisticsForPREDA <- StatisticsForPREDAFromdataframe(StatisticsForPREDA_dataframe=CNdataFile,ids_column=1,statistic_columns=2:(ncol(cn$data)+1), testedTail="both")

CNGenomicsAnnotations <- new("GenomicAnnotations", 
        ids = as.character(rownames(cn$data)), 
        chr = as.numeric(cn$info$chr), 
        start = as.numeric(cn$info$start), 
        end = as.numeric(cn$info$end),
        strand = as.numeric(rep(1,nrow(cn$data))), 
        chromosomesNumbers = chromosomes,
        chromosomesLabels = as.character(chromosomes))

CNGenomicsAnnotationsForPREDA <- GenomicAnnotations2GenomicAnnotationsForPREDA(CNGenomicsAnnotations, reference_position_type="median")

SODEGIRCNDataForPREDA <- MergeStatisticAnnotations2DataForPREDA(CNStatisticsForPREDA, CNGenomicsAnnotationsForPREDA, sortAndCleanNA = TRUE, quiet = FALSE, MedianCenter = TRUE)

# copy number analysis
SODEGIRCNanalysisResults <- PREDA_main(SODEGIRCNDataForPREDA, outputGenomicAnnotationsForPREDA = SODEGIRGEDataForPREDA, smoothMethod = smoothMethod, nperms = nperm)

SODEGIR_CN_GAIN <- PREDAResults2GenomicRegions(SODEGIRCNanalysisResults,
    qval.threshold = cn.qval.threshold, smoothStatistic.tail = "upper",
    smoothStatistic.threshold = cn.smoothStatistic.threshold.gain)
SODEGIR_CN_LOSS <- PREDAResults2GenomicRegions(SODEGIRCNanalysisResults,
    qval.threshold = cn.qval.threshold, smoothStatistic.tail = "lower",
    smoothStatistic.threshold = cn.smoothStatistic.threshold.loss)

# Check: right order for both data types?
#if(all(analysesNames(SODEGIRCNanalysisResults) == analysesNames(SODEGIRGEanalysisResults)) == FALSE){stop("Samples must be in same order for both gene expression and copy number data.")}

# joint SODEGIR analysis 
# analsze overlaps
SODEGIR_AMPLIFIED <- GenomicRegionsFindOverlap(SODEGIR_GE_UP, SODEGIR_CN_GAIN)
SODEGIR_DELETED <- GenomicRegionsFindOverlap(SODEGIR_GE_DOWN, SODEGIR_CN_LOSS)

names(SODEGIR_AMPLIFIED) <- names(SODEGIR_GE_UP)
names(SODEGIR_DELETED) <- names(SODEGIR_GE_DOWN)

# Calculate dataset signatures
SDGsignature_amplified <- computeDatasetSignature(SODEGIRGEDataForPREDA, genomicRegionsList = SODEGIR_AMPLIFIED)
SDGsignature_deleted <- computeDatasetSignature(SODEGIRGEDataForPREDA, genomicRegionsList = SODEGIR_DELETED)

# Compare found significant regions with
regions <- rbind(GenomicRegions2dataframe(SDGsignature_amplified[[1]]),GenomicRegions2dataframe(SDGsignature_deleted[[1]]))

found_index <- numeric()
for(j in 1:nrow(regions)){
    neu <- which(ge$info$chr == regions$chr[j] & ge$info$loc >= regions$start[j] & ge$info$loc <= regions$end[j])
    found_index <- c(found_index,neu)
}
all_genes <- rownames(ge$data)

found_genes <- all_genes[found_index]
not_found_genes <- setdiff(all_genes, found_genes)
false_negative_index <- which(all_genes %in% intersect(not_found_genes,cancerGenes))
false_negative_genes <- all_genes[false_negative_index]

not_cancerGenes <- setdiff(all_genes,cancerGenes)
false_positive_index <- which(all_genes %in% intersect(found_genes,not_cancerGenes))
false_positive_genes <- all_genes[false_positive_index]

true_positive_genes <- setdiff(found_genes,false_positive_genes)
true_negative_genes <- setdiff(not_found_genes,false_negative_genes)

ordg_best_case <- as.character(c(true_positive_genes,false_positive_genes,false_negative_genes,true_negative_genes))
ordg_worst_case <- as.character(c(false_positive_genes,true_positive_genes,true_negative_genes,false_negative_genes))

list(best_case_order=ordg_best_case, worst_case_order=ordg_worst_case)
}
