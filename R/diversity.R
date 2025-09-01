###############################################################################
# Diversity Measures Analysis 
###############################################################################
#
diversity.getDiversityMeasures <- function() {
    c("maxBarcodeFreq",
      "barcodeHet",
      "meanVariantHet",
      "meanHetCallsProp",
      "meanDistance",
      "medianDistance")
}
#
diversity.getDatasetForMeasure <- function(measure) {
    if (measure %in% c("maxBarcodeFreq", "barcodeHet")) { 
        datasetName <- "imputed"
    } else if (measure %in% c("meanVariantHet", "meanHetCallsProp")) { 
        datasetName <- "unfiltered"
    } else if (measure %in% c("meanDistance", "medianDistance")) { 
        datasetName <- "filtered"
    }
    datasetName
}
#
# Diversity measure estimates from genetic barcodes.
# Note that this is executed over imputed barcodes, and therefore there is not missingness or het genotypes in the barcodes.
#
diversity.estimateDiversityMeasures <- function (ctx, sampleNames, measureNames) {
    result <- c()
    for (mIdx in 1:length(measureNames)) {
        measureName <- measureNames[mIdx]			;print(paste("diversity.estimateDiversityMeasures", measureName))
        #
        # maxBarcodeFreq - this is the proportion of samples thaty carry the most common barcode (imputed)
        #
        if (measureName == "maxBarcodeFreq") {
	    value <- diversity.estimateMaxBarcodeFreq (ctx, sampleNames)
        } else if (measureName == "barcodeHet") {
	    value <- diversity.estimateBarcodeHet (ctx, sampleNames)
        } else if (measureName == "meanVariantHet") {
	    value <- diversity.estimateMeanVariantHet (ctx, sampleNames)
        } else if (measureName == "meanHetCallsProp") {
	    value <- diversity.estimateMeanHetCallsProp (ctx, sampleNames)
        } else if (measureName %in% c("meanDistance", "medianDistance")) {
	    value <- diversity.estimateAverageDistance (ctx, sampleNames, measureName)
        } else {
            stop(paste("Invalid diversity measure:", measureNames))
        }
        result <- c(result, value)
    }
    result
}
#
#
# maxBarcodeFreq - this is the proportion of samples thaty carry the most common barcode (imputed)
#
diversity.estimateMaxBarcodeFreq <- function (ctx, sampleNames) {
    allBarcodes <- ctx$rootCtx$imputed$barcodes
    barcodes <- allBarcodes[sampleNames]
    sampleCount <- length (barcodes)
    maxBcCounts <- max(table(barcodes))
    value <- maxBcCounts / sampleCount
    value
}
#
# barcodeHet - the probability that two samples picked at random will have different imputed haplotypes 
#            (i.e. the heterozygosity in this sample set, using imputed barcodes as the alleles).
#
diversity.estimateBarcodeHet <- function (ctx, sampleNames) {
    allBarcodes <- ctx$rootCtx$imputed$barcodes
    barcodes <- allBarcodes[sampleNames]
    sampleCount <- length (barcodes)
    probs <- table(barcodes) / sampleCount
    value <- (sampleCount / (sampleCount-1)) * (1 - sum(probs * probs))
    value
}
#
# meanVariantHet - the mean of the heterozygosity at each of the barcoding loci (using non-imputed data).
#              At a barcoding locus, heterozygosity is the probability that two samples picked at random will have different nucleotide alleles. 
#              Samples with missing genotype are ignored, while het samples are assigned "half" to each allele.
#
diversity.estimateMeanVariantHet <- function (ctx, sampleNames) {
    genoData <- ctx$rootCtx$unfiltered$barcodeGenoData
    varCount <- genoData$columnCount					#;print(paste("varCount", varCount))
    hets <- rep(NA, varCount)
    for (varIdx in 1:varCount) {					#;print(paste("varIdx", varIdx))
        alleleCounts <- list()
        snpGenos <- genoData$columnGenoData[[varIdx]]$sampleGenotypes
        snpAlleles <- genoData$columnGenoData[[varIdx]]$sampleAlleles
        for (sampleIdx in 1:length(sampleNames)) {				#;print(paste("sampleIdx", sampleIdx))
            sampleName <- sampleNames[sampleIdx]			#;print(paste("sampleName", sampleName))
            geno <- snpGenos[sampleName] 				#;print(paste("geno", geno))
            if (geno != GENO.MISS) {
                alleleProps <- snpAlleles[[sampleName]]			#;print(paste("alleleProps", alleleProps))
                alleles <- names(alleleProps)				#;print(paste("alleles", alleles))
                for (aIdx in 1:length(alleles)) {			#;print(paste("aIdx", aIdx))
                    a <- alleles[aIdx]					#;print(paste("a", a))
                    p <- alleleProps[[aIdx]]				#;print(paste("p", p))
                    if (is.null(alleleCounts[[a]])) {
                        alleleCounts[[a]] <- p
                    } else {
                        alleleCounts[[a]] <- alleleCounts[[a]] + p
                    }							#;print(paste("alleleCounts[[a]]", alleleCounts[[a]]))
                }
            }
        }
        aCounts <- unlist(alleleCounts)					#;print(paste("aCounts", aCounts))
        totalCount <- sum(aCounts)					#;print(paste("totalCount", totalCount))
        if (totalCount > 1) {
            aProps <- aCounts / totalCount				#;print(paste("aProps", aProps))
            aSq <- aProps * aProps
            hets[varIdx] <- (totalCount / (totalCount - 1)) * (1 - sum(aSq))	#;print(paste("hets[varIdx]", hets[varIdx]))
        }
    }
    value <- mean(hets, na.rm=TRUE)					#;print(paste("value", value))
    value
}
#
# meanHetCallsProp - the mean of the proportion of heterozygous calls at each of the barcoding loci (using non-imputed data).
#                  Samples with missing genotype are ignored.
#
diversity.estimateMeanHetCallsProp <- function (ctx, sampleNames) {
    genoData <- ctx$rootCtx$unfiltered$barcodeGenoData
    snpCount <- genoData$columnCount
    hetProps <- rep(NA, snpCount)
    for (snpIdx in 1:snpCount) {
        snpGenos <- genoData$columnGenoData[[snpIdx]]$sampleGenotypes
        genos <- snpGenos [sampleNames]
        homCount <- length(which(genos == GENO.HOM))
        nonMissingCount <- length(genos) - length(which(genos == GENO.MISS))
        if (nonMissingCount > 0) {
            hetProps[snpIdx] <- 1 - (homCount / nonMissingCount)
        }
    }
    value <- mean(hetProps, na.rm=TRUE)
    value
}
#
# medianDistance, meanDistance  - the median or mean of the pairwise genetic distance between the samples 
#
diversity.estimateAverageDistance <- function (ctx, sampleNames, measureName) {
    distData <- ctx$rootCtx$filtered$distData
    mat <- as.matrix(distData[sampleNames,sampleNames])
    mat[lower.tri(mat,diag=TRUE)] <- NA					#;print(mat)
    if (measureName == "medianDistance") {
        value <- stats::median(mat, na.rm=TRUE)
    } else if (measureName == "meanDistance") {
        value <- mean(mat, na.rm=TRUE)
    } 
    value
}
