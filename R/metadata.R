###############################################################################
# Metadata Column Statistics - General
################################################################################
#
# Returns an allele sample count summary string for the columns specified (one line per column).
# The value sample counts are semicolon-separated.
#
meta.getValueCounts <- function (ctx, sampleMeta, countableFeatureNames, params=NULL) {	#; print("meta.getValueCounts"); print(countableFeatureNames)
    if (is.null(countableFeatureNames) || (length(countableFeatureNames) == 0)) {
        return (c());
    }
    config <- context.getConfig(ctx)
    countableFeatures <- config$countableFeatures				#; print(head(countableFeatures))
    countableFeatures <- countableFeatures[countableFeatureNames,]		#; print(rownames(countableFeatures))
    colNames  <- setup.getFeatureColumns(countableFeatures)			#; print(colNames)
    result <- c()
    for (mIdx in 1:length(colNames)) {
        colName <- colNames[mIdx]						#; print(colName)
        counts <- meta.getColumnValueCounts (ctx, sampleMeta, colName)
        nVals <- names(counts)
        valPairs <- paste(nVals, counts, sep=":")
        valPairStr <- paste(valPairs, collapse="; ")
        result <- c(result, valPairStr)						#; print(valPairStr)
    }
    result <- as.character(result)
    names(result) <- countableFeatureNames
    result
}
#
# Returns a tabulated array of alleles sample counts in a column, sorted by decreasing sample count
# The names of the elements are the allele labels. Missing and <NA> are discarded.
#
meta.getColumnValueCounts <- function (ctx, sampleMeta, columnName) {	#; print("meta.getColumnValueCounts"); print(columnName)
    sampleNames <- rownames(sampleMeta)					#; print(length(sampleNames))
    columnGenos <- ctx$rootCtx$alleleGenoData$columnGenoData[[columnName]]
    sampleGenos <- columnGenos$sampleGenotypes[sampleNames]
    
    vals <- c()
    for (sIdx in 1:length(sampleNames)) {
        sampleGeno <- sampleGenos[sIdx]
        if (sampleGeno==GENO.HOM) {
            sName <- sampleNames[sIdx]
            sAlleleProps <- columnGenos$sampleAlleles[[sName]]
            sVals <- names(sAlleleProps)
            vals <- c(vals, sVals[1])
        }
    }
    counts <- table(vals)						#; print(counts)
    counts <- sort(counts, decreasing=TRUE)
    counts
}
#
###############################################################################
# Metadata Statistics - Resistance Predictions
################################################################################
#
#
#
meta.getResistancePrevalence <- function (ctx, sampleMeta, drugNames, params=NULL) {	#; print(nrow(sampleMeta))
    if (is.null(drugNames) || (length(drugNames) == 0)) {
        return (c());
    }
    config <- context.getConfig(ctx)
    aggregateCountMin <- 1
    if (!is.null(params)) {
        aggregateCountMin <- param.getParam ("map.aggregateCountMin", params)
    }  							#; print(aggregateCountMin)
    result <- c()
    for (mIdx in 1:length(drugNames)) {
        # Get the feature name, and the corresponding column in the GRC data
        drugName <- drugNames[mIdx]			#; print(drugName)
        feat <- config$drugPredictionFeatures
        drugCol <- feat[drugName, "ColumnName"]		#; print(drugCol)
        phenos <- sampleMeta[,drugCol]			#; print(phenos)
        rCount <- length(which(phenos=="Resistant"))	#; print(rCount)
        sCount <- length(which(phenos=="Sensitive"))	#; print(sCount)
        totalCount <- rCount + sCount
        if (totalCount >= aggregateCountMin) {
            preval <- rCount / totalCount
        } else {
            preval <- NA
        }
        result <- c(result, preval)			#; print(preval)
    }
    result <- as.numeric(result)			#; print(result)
    names(result) <- drugNames				#; print(result)
    result
}
#
###############################################################################
# Metadata Statistics - V2
################################################################################
#
# Mutation Prevalence
#
meta.getMutationPrevalence <- function (ctx, sampleMeta, mutationNames, params) {	#; print("meta.getMutationPrevalence")#; print(mutationNames)
    config <- context.getConfig(ctx)
    aggregateCountMin <- 1
    if (!is.null(params)) {
        aggregateCountMin <- param.getParam ("map.aggregateCountMin", params)
    }  								#; print(aggregateCountMin)
    prevalenceData <- c()
    for (mIdx in 1:length(mutationNames)) {
        mutName <- mutationNames[mIdx]				#; print(mutName)
        feat <- config$drugMutationFeatures
	columnName <- feat[mutName, "ColumnName"]		#; print(columnName)
        mutAllele <- feat[mutName, "Alternative"]		#; print(mutAllele)
        #
        #mutData <- meta.getPositionData (config, mutName)
        #if (is.null(mutData)) {
        #    stop(paste("Invalid mutation specified:",mutName))
        #}
        #
        sampleNames <- rownames(sampleMeta)			#; print(length(sampleNames))
        columnGenos <- ctx$rootCtx$alleleGenoData$columnGenoData[[columnName]]
        sampleGenos <- columnGenos$sampleGenotypes[sampleNames]
        prop <- 0.0; count <- 0
        for (sIdx in 1:length(sampleNames)) {
            sampleGeno <- sampleGenos[sIdx]
            if ((sampleGeno==GENO.MISS)||(sampleGeno==GENO.HET_NO_PROPS)) {
                next
            }
            sName <- sampleNames[sIdx]
            sAlleleProps <- columnGenos$sampleAlleles[[sName]]
            propMut <- sAlleleProps[[mutAllele]]
            if (!is.null(propMut)) {
                prop <- prop + propMut
            }
            count <- count + 1
        }
        prevalence <- NA
        if (count >= aggregateCountMin) {
            prevalence <- prop / count
        }
        prevalenceData <- c(prevalenceData, prevalence)		#; print(prevalence)
    }
    prevalenceData <- as.numeric(prevalenceData)
    names(prevalenceData) <- paste(mutationNames)
    prevalenceData
}
#
###############################################################################
# Metadata processing
################################################################################
#
# Select sample metadata as specified in the config
#
meta.select <- function(sampleMeta, selectDefs) {
    for (sIdx in 1:length(selectDefs)) {
        selDef <- selectDefs[[sIdx]]
        selField <- selDef$field
        selValues <- selDef$values
        #print (paste("Select field: ",selField))
        #print (paste("Select values: ",selValues))
        metaValues <- sampleMeta[,selField]
        sampleMeta <- sampleMeta[which (metaValues %in% selValues),]
        #print (paste("Selected samples: ",nrow(sampleMeta)))
    }
    sampleMeta
}
#
# Select sample metadata within a given time period (start and end date inclusive)
#
meta.defaultDateFormat <- "%d/%m/%Y"
meta.naDate <- "01/01/1000"		# Date value to use instead of NA or "-"
#
meta.filterByDate <- function(sampleMeta, startDate, endDate, format=meta.defaultDateFormat) {
    metaDates <- meta.getSampleDates (sampleMeta, format)	#; print(head(metaDates))
    selectIdx <- !is.na(metaDates)
    if (!is.null(startDate)) {
        selectIdx <- selectIdx & (metaDates>=startDate)
    }								#; print(selectIdx[1:10])
    if (!is.null(endDate)) {
        selectIdx <- selectIdx & (metaDates<=endDate)
    }								#; print(selectIdx[1:10])
    filterMeta <- NULL						#; print(length(which(selectIdx==TRUE)))
    if (length(which(selectIdx==TRUE)) > 0) {
         filterMeta <- sampleMeta[selectIdx,]
    }								#; print(nrow(filterMeta))
    filterMeta
}
#
#
#
meta.getSampleDates <- function (sampleMeta, format=meta.defaultDateFormat) {
    .naDate <- 
    sampleCount <- nrow(sampleMeta)
    if ("CollectionDate" %in% colnames(sampleMeta)) {
        metaDatesIn <- as.character(sampleMeta$CollectionDate)
    } else {
        #message("Your GRC data does not contain the CollectionDate Field. Using the Year field instead- see the documentation.")
        metaDatesIn <- rep("-", sampleCount)
    }									#; print(metaDatesIn[1:500])
    #
    # If date is missing, use the year which shoul be there
    #
    missingDateIndexes <- which (metaDatesIn == "-")			#; print(missingDateIndexes)
    if (length(missingDateIndexes) > 0) {
        years <- suppressWarnings(as.integer(sampleMeta$Year))		#; print(years[1:500])
        yearDates <- paste0("31/12/", years)				#; print(yearDates[1:500])
        for (mdIdx in missingDateIndexes) {
            if (is.na(years[mdIdx])) {
                #
                # Silly workaround. Placing a "-" can cause as.Date() errors if it is in the first element of the vector!!!
                # So we put a valid but improbably date instead. Not ideal, but avoids errors
                #
                #metaDatesIn[mdIdx] <- "-" 
                metaDatesIn[mdIdx] <- meta.naDate 
            } else {
                metaDatesIn[mdIdx] <- yearDates[mdIdx]
            }
        }								#; print(metaDatesIn[1:500])
    }
    metaDates <- as.Date(metaDatesIn, tryFormats=format, optional=TRUE)	#; print(metaDates[1:500])
    metaDates[which(metaDates==meta.naDate)] <- NA
    metaDates
}
