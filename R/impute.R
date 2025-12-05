###################################################################
# Imputed Data, with filled-in missingness in filtered barcodes 
###################################################################
#
impute.DEFAULT_MOST_SIMILAR_COUNT <- 100
#
#
#
impute.selectImputableSamples <- function (ctx, barcodeGenoData) {
    config        <- context.getConfig(ctx)
    maxImputeProp <- config$maxImputedProportion
    print("Filtering barcode genotypes by imputing proportion required")
    colCount      <- barcodeGenoData$columnCount
    missingCounts <- barcodeGenoData$sampleMissingCounts
    hetCounts     <- barcodeGenoData$sampleHetCounts
    imputeProps   <- (missingCounts + hetCounts) / colCount
    samples       <- barcodeGenoData$samples
    selSamples    <- samples[which(imputeProps <= maxImputeProp)]
    selSamples
}
#
#
#
impute.getMostCommonAlleleInColumn <- function (columnGenos) {
    sGenotypes <- columnGenos$sampleGenotypes				#; print(length(sGenotypes))
    sAlleles   <- columnGenos$sampleAlleles				#; print(length(sAlleles))
    sampleCount <- length(sGenotypes)
    homAlleles <- rep("-", sampleCount)
    #
    for (sIdx in 1:sampleCount) {
        if (sGenotypes[sIdx] == GENO.HOM) {
            homAlleles[sIdx] = names(sAlleles[[sIdx]])
        }
     }
     homAlleles <- homAlleles[which(homAlleles != "-")]
     homAlleleCounts <- sort(table(homAlleles), decreasing=TRUE)	#; print(homAlleleCounts)
     mostCommonAllele <- names(homAlleleCounts)[1]
     mostCommonAllele
}
#
#
#
impute.imputeBarcodeGenos <- function (ctx, barcodeGenoData) {
    sampleNames <- barcodeGenoData$samples
    sampleCount <- length(sampleNames)						#; print(sampleCount)
    msCount <- min((sampleCount-1), impute.DEFAULT_MOST_SIMILAR_COUNT)		#; print(msCount)		
    # 
    # Now get the data for the 100 most similar samples to each sample in the matrix.
    # Arrange the distance matrix so that the sample name sequence is that of the genotype data,
    # or else the sample indexes or the "most similar sample" tables will be scrambled
    #
    distData <- context.getDistanceMatrix (ctx, sampleSetName=NULL, useImputation=FALSE)	#; print(dim(distData))
    distData <- distData[sampleNames,sampleNames]				#; print(dim(distData))
    msData   <- impute.findMostSimilarSamples(distData, msCount)		#; print(dim(msData))
    msIndexData    <- msData$indexes						#; print(dim(msIndexData))
    msDistanceData <- msData$distances						#; print(dim(msIndexData))
    #
    # Get some vectors for the allele analyses
    #
    msAlleles   <- character(msCount)
    msDistances <- numeric(msCount)
    #
    #
    #
    newColumnGenoData <- list()
    #
    columnGenoData <- barcodeGenoData$columnGenoData
    columnNames <- names(columnGenoData)					#; print(head(columnNames))
    columnCount <- length(columnNames)						#; print(columnCount)
    for (vIdx in 1:columnCount) {						#; print(vIdx)
        mostCommonAllele <- NULL
        columnGenos <- columnGenoData[[vIdx]]					#; print("vIdx"); print(vIdx)
        sGenotypes <- columnGenos$sampleGenotypes				#; print(sGenotypes)
        sAlleles   <- columnGenos$sampleAlleles					#; print(length(sAlleles))
        newSampleAlleles <-list()
        for (sIdx in 1:sampleCount) {						#; print(sIdx)
            #
            # If geno is a HOM, no need to impute, transfer the current allele
            #
            if (sGenotypes[sIdx] == GENO.HOM) {
                newSampleAlleles[[sIdx]] <- sAlleles[[sIdx]]
                next
            }
            #
            # We need to impute, by looking at the alleles in the 100 nearest samples
            #
            msSampleIndexes <- msIndexData[,sIdx]				#; print(msSampleIndexes)
            msSampleDistances <- msDistanceData[,sIdx]
            for (msIdx in 1:msCount) {
                msSampleIndex <- msSampleIndexes[msIdx]				#; print(paste(msIdx, "-", msSampleIndex))
                if (sGenotypes[msSampleIndex] == GENO.HOM) {
                    msAlleles[msIdx] <- names(sAlleles[[msSampleIndex]])
                    msDistances[msIdx] <- msSampleDistances[msIdx]
                } else {
                    msAlleles[msIdx] <- "-"
                }
            }									#; print(msAlleles)
            #
            nonMissIdx <- which(msAlleles != "-")
            if (length(nonMissIdx) == 0) {
                #
                # Special case: all 100 nearest samples are missing genotypes.
                # This could mean there's a deletion or something, but here we 
                # deal with it simply by imputing the commonest allele.
                #
                if (is.null(mostCommonAllele)) {				#; print(paste(vIdx, sIdx))
                    mostCommonAllele <- impute.getMostCommonAlleleInColumn (columnGenos)
                }
                allele <- mostCommonAllele					#; print(allele)
            } else {
                #
                # Work out the allele scores, and pick the one with the lowest
                # To score the alleles, we use a distance score, which assigns to each allele a value between 0 and 1, 
                # depending on the genetic distance of their sample from the sample being imputed. 
                # The score is then squared to exaggerate the contribution of samples that are very close.
                #
                msd <- msDistances[nonMissIdx]					#; print(msd)
                msa <- msAlleles[nonMissIdx]					#; print(msa)
                names(msd) <- msa						#; print(msd)
                msd <- sort(msd)						#; print(msd)
                d1 <- msd[1]							#; print(d1)
                dk <- msd[length(msd)]						#; print(dk)
                range <- (dk - d1)						#; print(range)
                if (range > 0.0) {
                    weights <- (dk - msd) / range
                    weights <- weights * weights; 	# Note: the weights are squared!
                } else {
                     weights <- rep(1, length(msd))
                }
                names(weights) <- names(msd)					#; print(weights)
                table <- tapply(weights, names(weights), FUN=sum)		#; print(table)
                sorted <- sort(table, decreasing=TRUE)				#; print(sorted)
                allele <- names(sorted)[1]            				#; print(allele)
            }

            #
            # Change the alleles as appropriate
            #
            aList <- list(1); names(aList) <- allele
            newSampleAlleles[[sIdx]] <- aList					#; print(newSampleAlleles[[sIdx]])
            #print("end imputeBarcodeGenos")
        }
        newSampleGenotypes <- rep(GENO.HOM, sampleCount)
        newColGenotypes <- list(samples=sampleNames,
                                sampleGenotypes=newSampleGenotypes,
                                sampleAlleles=newSampleAlleles,
                                totalCount=sampleCount,
                                missingCount=0,
                                homCount=sampleCount,
                                hetCount=0)
        colName <- columnNames[vIdx]
        newColumnGenoData[[colName]] <- newColGenotypes
    }
    #
    # Estimate sample missingness and heterozyous call proportion
    #
    sProp <- geno_estimateSampleProperties (sampleNames, columnNames, newColumnGenoData)
    cProp <- geno_estimateColumnProperties (columnNames, newColumnGenoData)
    #
    newGenotypeData <- list(samples=sampleNames, 
                            columns=columnNames,
                            columnGenoData=newColumnGenoData, 
                            #
                            columnCount=sProp$columnCount, 
                            sampleMissingCounts=sProp$missingCounts, 
                            sampleHomCounts=sProp$homCounts,
                            sampleHetCounts=sProp$hetCounts,
                            #
                            sampleCount=cProp$sampleCount, 
                            columnMissingCounts=cProp$missingCounts, 
                            columnHomCounts=cProp$homCounts, 
                            columnHetCounts=cProp$hetCounts
                            )
    newGenotypeData
}
#
#
#
impute.findMostSimilarSamples <- function (distData, mostSimilarCount) {
    #
    # Get the indexes of all the samples
    #
    sampleNames <- colnames(distData)
    sampleCount <- length(sampleNames)
    sampleIndexes <- 1:sampleCount
    names(sampleIndexes) <- sampleNames
    #
    # Create data frames for the results
    #
    msIdxData  <- data.frame(matrix(nrow=mostSimilarCount,ncol=0))
    msDistData <- data.frame(matrix(nrow=mostSimilarCount,ncol=0))
    #
    # Create two columns for each sample, and attache them 
    #
    for (sIdx in 1:sampleCount) {				#; isHere <- (sIdx==398); print("findMostSimilarSamples")
        sampleName <- sampleNames[sIdx]				#; print(paste(sIdx,sampleNames[sIdx]))
        #
        # Get the distances of each sample vs the test sample, removing the test sample itself
        #
        sDist <- distData[,sIdx]
        names(sDist) <- sampleNames				#; print(length(sDist))
        sDist <- sDist[which(sampleNames != sampleName)]	#; print(length(sDist)); print(sDist)
        #
        # Pick the 100 samples with the least distance
        #
        sDist <- sort(sDist)
        msDist <- sDist[1:mostSimilarCount]			#; print(msDist)
        msSamples <- names(msDist)				#; print(msSamples)
        msIdxData[,sampleName] <- sampleIndexes[msSamples]	#; print(msIdxData[,sampleName])
        msDistData[,sampleName] <- msDist			#; print(msDistData[,sampleName]); print("end findMostSimilarSamples")
    }
    list(indexes=msIdxData, distances=msDistData)
}
#
#
#
impute.buildBarcodeGenotypeTable <- function (barcodeGenoData) {
    sampleNames <- barcodeGenoData$samples
    sampleCount <- length(sampleNames)			#; print(sampleCount)
    columnGenoData <- barcodeGenoData$columnGenoData
    columnNames <- names(columnGenoData)		#; print(head(columnNames))
    columnCount <- length(columnNames)			#; print(columnCount)
    barcodeGenoTable <- data.frame(matrix(nrow=sampleCount, ncol=0))
    sGenos <- character(sampleCount)
    for (vIdx in 1:columnCount) {
        columnGenos <- columnGenoData[[vIdx]]
        sampleAlleles <- columnGenos$sampleAlleles	#; print(length(sampleAlleles))
        for (sIdx in 1:sampleCount) {
            sGenos[sIdx] = names(sampleAlleles[[sIdx]])
        }
        #alleleCounts <- sort(table(sGenos), decreasing=TRUE)	#; print(alleleCounts)
        #alleles <- names(alleleCounts)
        #alleleIndexes <- 1:length(alleles)
        #names(alleleIndexes) <- alleles
        #sGenoIndexes <- alleleIndexes[sGenos]
        #barcodeGenoTable <- cbind(barcodeGenoTable, sGenoIndexes)
        barcodeGenoTable <- cbind(barcodeGenoTable, sGenos)
    }
    rownames(barcodeGenoTable) <- sampleNames
    colnames(barcodeGenoTable) <- columnNames
    barcodeGenoTable
}
#
#
#
impute.buildBarcodeStrings <- function (barcodeGenoTable) {
    sampleNames <- rownames(barcodeGenoTable)
    sampleCount <- length(sampleNames)
    barcodes <- character(sampleCount)
    for (sIdx in 1:sampleCount) {
        barcodes[sIdx] <- paste(barcodeGenoTable[sIdx,], collapse="-")
    }
    names(barcodes) <- sampleNames
    barcodes
}
#
#
#
impute.getMostCommonAlleleInColumn <- function (columnGenos) {
    sGenotypes <- columnGenos$sampleGenotypes				#; print(length(sGenotypes))
    sAlleles   <- columnGenos$sampleAlleles				#; print(length(sAlleles))
    sampleCount <- length(sGenotypes)
    homAlleles <- rep("-", sampleCount)
    #
    for (sIdx in 1:sampleCount) {
        if (sGenotypes[sIdx] == GENO.HOM) {
            homAlleles[sIdx] = names(sAlleles[[sIdx]])
        }
     }
     homAlleles <- homAlleles[which(homAlleles != "-")]
     homAlleleCounts <- sort(table(homAlleles), decreasing=TRUE)	#; print(homAlleleCounts)
     mostCommonAllele <- names(homAlleleCounts)[1]
     mostCommonAllele
}
