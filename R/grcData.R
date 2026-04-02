###############################################################################
# Read a GRC data file from an Excel sheet
################################################################################
#
# Read sample data file
#
grcData.load <- function (grcDataFile, species, version, grcDataSheet) {
    if (!file.exists(grcDataFile)) {
        stop(paste("GRC data file", grcDataFile, "does not exist."))
    }
    
    if (hasExcelSheet (grcDataFile, "Properties")) {
        #
        # Read the GRC file properties, since this spreadsheet has one
        #
        if (!is.null(grcDataSheet)) {
            rlang::warn("Ignoring argument \"sheet\" that is only valid for versions 1.x.")
        }
        grcProps <- grcData.getGrcProperties (grcDataFile)
        species <- grcProps$species
        version <- grcProps$version
        provider <- grcProps$provider
        providerVersion <- grcProps$providerVersion
        #
        # Read the GRC file features
        #
        grcFeatures <- readExcelData(grcDataFile, "Features")
        grcFeatures <- grcData.checkColumns (grcFeatures, c("FeatureName","ColumnName","Class","DataType","WT","Alternative"), "Features", trim=TRUE)
        rownames(grcFeatures) <- grcFeatures$FeatureName
        grcDataSheet <- "Data"

    } else {
        #
        # No GRC file properties, so assume it's a V1.x GRC, the version has to be specified
        #
        if (is.null(version)) {
            stop(paste("You must specify a version number of the GRC data file."))
        }
        verParts <- unlist(strsplit(version,"\\."))
        majorVersion <- verParts[1]
        minorVersion <- verParts[2]
        if (majorVersion != 1) {
            stop(paste("GRC file without a Properties sheet can only be V1.x"))
        }
        #
        # The way in which this version checking is written is a bit odd, for historical reasons.
        # In V1, different species had different GRC format versions. Starting with V2, we have a GRC format version, used by all species.
        # Currently we still want to support V1.4 of Pf and 1.0 of Pv. This may change later.
        #
        if (is.null(species)) {
            stop("Species was not specified for the GRC data file")
        }
        if (species == "Pf") {
            if (version != "1.4") {
                stop(paste("Invalid species/version combination:",species, "/", version))
            }
            # The GRC file features are hardcoded
            grcFeatures <- grcData.v1.Pf.getFeatures()
            
        } else if (species == "Pv") {
            if (version != "1.0") {
                stop(paste("Invalid species/version combination:",species, "/", version))
            }
            # The GRC file features are hardcoded
            grcFeatures <- grcData.v1.Pv.getFeatures()
            
        } else {
            stop(paste("Invalid species:",species))
        }
        #
        # Assume only GenRe created V1.x GRCs.
        #
        provider <- "GenRe-Mekong"
        providerVersion <- version
    }
    grc <- NULL
    #
    # Read the GRC Sample Data
    # Make sure all the feature columns are present
    #
    grcSampleData <- readExcelData(grcDataFile, grcDataSheet)
    grcSampleData <- grcData.checkColumns (grcSampleData, c("SampleId",grcFeatures$ColumnName), grcDataSheet)
    rownames(grcSampleData) <- grcSampleData$SampleId	
    #
    # Ensure we have a column named CollectionDate that has a valid cate for each sample, determined from Year and Month if necessary
    #
    grcSampleData <- grcData.determineSampleDates(grcSampleData)
    #
    # Process the features to create a configuation that works for this GRC
    #
    grcConfig <- grcData.buildGrcConfig (species, version, grcFeatures)
    #
    grc <- list(data=grcSampleData, config=grcConfig, version=version, provider=provider, providerVersion=providerVersion)
    grc
}
#
#
#
grcData.getGrcProperties <- function (grcDataFile) {
    propsData <- readExcelData(grcDataFile, "Properties")
    propsData <- grcData.checkColumns (propsData, c("Name","Value"), "Properties", trim=TRUE)
    props        <- propsData$Value
    names(props) <- propsData$Name
    #
    fileFormat <- grcData.getGrcPropertyValue ("Format", props, "GRC")
    majorVersion <- grcData.getGrcPropertyValue ("MajorVersion", props)
    minorVersion <- grcData.getGrcPropertyValue ("MinorVersion", props)
    version <- paste(majorVersion, minorVersion, sep=".")
    if (as.integer(majorVersion) > 2) {
        stop(paste0("Unsupported GRC version: ", version))
    }
    if (as.integer(minorVersion) > 0) {
        stop(paste0("Unsupported GRC version: ", version))
    }
    provider <- grcData.getGrcPropertyValue ("Provider", props)
    providerVersion <- grcData.getGrcPropertyValue ("ProviderVersion", props)
    #
    supportedSpecies <- c(
        "Plasmodium falciparum"="Pf", 
        "Plasmodium vivax"="Pv"
    )
    speciesStr <- grcData.getGrcPropertyValue ("Species", props, names(supportedSpecies))
    species <- supportedSpecies[speciesStr]
    #
    grcProps <- list(species=species, version=version, provider=provider, providerVersion=providerVersion)
    grcProps
}
#
grcData.checkColumns <- function (data, columnNames, sheetName, trim=FALSE) {
    missing <- c()
    for (col in columnNames) {
        if (!(col %in% colnames(data))) {
            missing <- c(missing, col)
        } else {
            if (trim) {
                data[,col] <- trimws(data[,col])
            }
        }
    }
    if (length(missing) > 0) {
        missFields <- paste(missing, sep=",")
        stop(paste0("Missing column(s) in GRC sheet '",sheetName,"': ",missFields))
    }
    data
}
#
grcData.getGrcPropertyValue <- function (propName, props, expected=NULL) {
    if (!(propName %in% names(props))) {
        stop(paste("Missing GRC file property:", propName))
    }
    value <- props[propName]
    if (!is.null(expected)) {
        if (!(value %in% expected)) {
            stop(paste0("Invalid value spectified GRC file property '",propName,"': ", value))
        }
    }
    value
}
#
# #####################################################################################
# Sample CollectionDate processing.
# Ensure we have a column named DeterminedDate that has a valid date for each sample, determined from Year and Month if necessary.
# 
# GRC Sample Date rules:
# - For any sample, CollectionDate may be missing, either because it is missing for the sample, or because the CollectionDate column is missing (Public Release)
# - When present, CollectionDate must be a String in format "dd/mm/yyyy". This must be done in the Metadata file.
# - The Year value is mandatory; if the CollectionDate is provided, they must be consistent. This must be done in the Metadata file.
# - The Month value is optional; if the CollectionDate is provided, they must be consistent. This must be done in the Metadata file.
# - As for other columns, missing value is "-"
#
# Selection of samples by date intervals will be based on a metadata column called DeterminedDate, derived from Year, Month and CollectionDate
#
# #####################################################################################
#
grcData.determineSampleDates <- function(grcSampleData) {
    sampleNames <- grcSampleData$SampleId
    sampleCount <- nrow(grcSampleData)
    #
    # Get the input columns from the metadata
    #    
    if ("CollectionDate" %in% colnames(grcSampleData)) {
        collDate <- as.character(grcSampleData$CollectionDate)
    } else {
        collDate <- rep("-", sampleCount)
    }
    if ("Month" %in% colnames(grcSampleData)) {
        collMonth <- suppressWarnings(as.integer(grcSampleData$Month))
    } else {
        collMonth <- rep(NA_integer_, sampleCount)
    }
    collYear <- suppressWarnings(as.integer(grcSampleData$Year))
    badYearIdx <- which(is.na(collYear))
    if (length(badYearIdx) > 0) {
        badYearSampleIds <- sampleNames[badYearIdx]
        warning(paste("Invalid 'Year' value found in the following sample(s):", paste(badYearSampleIds, collapse=",")))
    }
    #
    # Go through each sample, determining the dates
    #
    detDates <- c()
    for (sIdx in 1:sampleCount) {
        if (collDate[sIdx] != "-") {
            # If the CollectionDate is missing, work one out from Year and Month
            year <- collYear[sIdx]
            month <- collMonth[sIdx]
            if (is.na(month)) {
                month <- 7
                day <- 1
            } else {
                day <- 15
            }
            cd <- lubridate::make_date(year=year, month=month, day=day)
        } else {
            # Is the CollectionDate date in the right format?
            cd <- as.Date(sIdx, "%d/%m/%Y")
            if (is.na(cd)) {
                warning(paste("Invalid 'CollectionDate' value found:", sampleNames[sIdx]))
            } else {
                if (lubridate::year(cd) != collYear[sIdx]) {
                    warning(paste("Inconsistency between 'CollectionDate' and 'Year' in sample", sampleNames[sIdx]))
                }
                if (lubridate::month(cd) != collMonth[sIdx]) {
                    warning(paste("Inconsistency between 'CollectionDate' and 'Month' in sample", sampleNames[sIdx]))
                }
            }
        }
        detDates <- c(detDates, cd)
    }
    grcSampleData$DeterminedDate <- detDates
    grcSampleData
}
#
# #####################################################################################
# GRC-specific setup configurations
# #####################################################################################
grcData.buildGrcConfig <- function (species, version, grcFeatures) {
    #
    config <- list (
        species = species,
        version = version
    )
    #
    # Drug Resistant Predictions
    #
    config$drugPredictionFeatures <- grcFeatures[which(grcFeatures$Class == "ResistancePrediction"),]
    #
    # Drug Resistant Positions
    #
    phenoFeatureClasses <- c("PhenotypeAssociatedLocus", "PhenotypeAssociatedVariant")
    config$drugLocusFeatures <- grcFeatures[which(grcFeatures$Class %in% phenoFeatureClasses),]	#; print(config$drugLocusFeatures)
    #
    # Drug Resistant Mutations
    #
    config$drugMutationFeatures <- grcFeatures[which(grcFeatures$Class == "PhenotypeAssociatedMutation"),]	#; print(config$drugMutationFeatures)
    #
    # Count Columns (i.e. columns used for Pie Charts)
    # excluding Amplifications
    #
    countColumns <- grcFeatures[which(grcFeatures$Class %in% c("MutationList","PhenotypeAssociatedLocus")),]
    #config$countableFeatures <- countColumns[which(countColumns$DataType != "Amplification"),]
    config$countableFeatures <- countColumns
    #
    # Amplifications Columns (also used for Pie Charts)
    #
    config$amplificationFeatures <- grcFeatures[which(grcFeatures$DataType == "Amplification"),]
    #
    # Barcoding Columns (also used for Pie Charts)
    #
    config$barcodingFeatures <- grcFeatures[which(grcFeatures$Class == "BarcodingLocus"),]
    #
    # Get the names of the allele-based columns to be genotyped, separating the barcoding SNPs
    #
    config$alleleColumns <- unique(c(config$drugLocusFeatures$ColumnName, 
                                     config$drugLocusFeatures$ColumnName, 
                                     config$drugMutationFeatures$ColumnName, 
                                     config$countableFeatures$ColumnName, 
                                     config$amplificationFeatures$ColumnName))
    config$barcodeColumns <- unique(config$barcodingFeatures$ColumnName)
    #
    # Copy the features used in clustering computations
    #
    config$cluster.stats.drugPredictionFeatures  <- config$drugPredictionFeatures
    config$cluster.stats.drugMutationFeatures    <- config$drugMutationFeatures
    config$cluster.stats.countableFeatures       <- config$countableFeatures
    #
    config
}
#
################################################################################
# Merge two GRC datasets
################################################################################
#
# Merge two GRC datasets
#
grcData.merge <- function (inGrc, newGrc, overwrite=FALSE, extendColumns=FALSE) {
    if ((inGrc$config$species != newGrc$config$species)) {
        stop("You can only merge GRCs for the same species")
    }
    if ((inGrc$provider != newGrc$provider) || (inGrc$providerVersion != newGrc$providerVersion)) {
        stop("Currently you can only merge GRCs from the same provider and with the same version number")
    }
    inData  <- inGrc$data
    newData <- newGrc$data

    inSamples  <- rownames(inData);	inFlds  <- colnames(inData)
    newSamples <- rownames(newData);	newFlds <- colnames(newData)
    #
    # Add missing columns in the new data if necessary
    #
    missNewFlds <- inFlds[which(!(inFlds %in% newFlds))]
    if (length(missNewFlds) > 0) {
        missNewData <- data.frame(matrix("<NA>", nrow=length(newSamples), ncol=length(missNewFlds)))
        colnames(missNewData) <- missNewFlds
        rownames(missNewData) <- newSamples
        newData <- cbind(newData, missNewData)
    }
    #
    # Remove extra columns if needed, or add columns introduced by the new dataset
    #
    extraNewFlds <- newFlds[which(!(newFlds %in% inFlds))]
    if (length(extraNewFlds) > 0) {
        if (extendColumns) {
            extraInData <- data.frame(matrix("<NA>", nrow=nrow(inSamples), ncol=length(extraNewFlds)))
            colnames(extraInData) <- extraNewFlds
            rownames(extraInData) <- inSamples
            inData <- cbind(inData, extraInData)
        } else {
            newData <- newData[,-extraNewFlds]
        }
    }
    #
    # Align the columns (they are the same in both datasets)
    #
    inFlds  <- colnames(inData)    
    newData <- newData[,inFlds]
    #
    # Remove overlap in samples before merging
    #
    overlap <- inSamples[which(inSamples %in% newSamples)]
    if (length(overlap) > 0) {
        if (overwrite) {
            inData <- inData[-overlap,]
        } else {
            newData <- newData[-overlap,]
        }
    }
    #
    # Finally, merge the rows
    #
    mergedData <- rbind(inData, newData)
    mergedGrc <- list(data=mergedData, config=inGrc$config, version=inGrc$version, 
                      provider=inGrc$provider, providerVersion=inGrc$providerVersion)
    mergedGrc
}
#
################################################################################
# GRC feature names listings
################################################################################
#
# 
#
grcData.getFeatureNames <- function (grc, featureType) {
    featureNames <- NULL
    config <- grc$config
    
    if (featureType == "drug") {
        featureNames <- config$drugPredictionFeatures$FeatureName
        
    } else if (featureType == "locus") {
        featureNames <- config$drugLocusFeatures$FeatureName
    
    } else if (featureType == "mutation") {
        featureNames <- config$drugMutationFeatures$FeatureName
    
    } else if (featureType == "amplification") {
        featureNames <- config$amplificationFeatures$FeatureName
    
    }
    featureNames
}
