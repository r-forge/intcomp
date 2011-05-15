process.ge <- function (ge, probespanGE) {

  # Convert gene expression data into mrnaSet object required by intCNGEan

  print("MRNA")

  # Ensure that sex chromosomes are numerical
  ychr <- ge$info$chr
  ychr <- gsub("X",23,ychr)
  ychr <- gsub("Y",24,ychr)

  df <- data.frame(Chromosome = as.numeric(as.character(ychr)),
	           Start = as.integer(ge$info$loc) - probespanGE,
		   End = as.integer(ge$info$loc) + probespanGE)

  featureData_ge <- new("AnnotatedDataFrame",
                        varMetadata = data.frame(labelDescription = c("Chromosomal position","Basepair position start","Basepair position end"),
                         row.names = c("Chromosome", "Start", "End")),dimLabels = c("featureNames", "featureColumns"), data = df)

  featureNames(featureData_ge) <- rownames(ge$data)
  mrnaSet <- new("ExpressionSet", exprs = as.matrix(ge$data), featureData = featureData_ge)

  mrnaSet

}