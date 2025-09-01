###############################################################################
# Map Aggregated Measure Analysis
###############################################################################
#
markerMap.executeMap <- function(map) {						#; print("markerMap.executeMap")
    mapMaster   <- map$master
    mapType     <- mapMaster$type						#; print(paste("mapType",mapType))
    measureName <- map$measureName						#; print(paste("measureName",measureName))
    interval    <- map$interval							#; print(interval)
    #
    datasetName <- map$datasetName						#; print(paste("datasetName",datasetName))
    sampleSet   <- mapMaster$sampleSet
    userCtx     <- mapMaster$userCtx
    params      <- mapMaster$params
    config      <- context.getConfig(userCtx)
    #
    sampleMeta <- context.getMeta (map$mapCtx, datasetName)			#; print(paste("markerMap.executeMap", nrow(sampleMeta)))
    if (nrow(sampleMeta)==0) {
        print(paste("No samples found - skipping interval", interval$name))
        return()
    }
    #
    baseMapInfo <- mapMaster$baseMapInfo
    #
    # Get the output folders
    #
    dataFolder <- mapMaster$dataFolder
    plotFolder <- mapMaster$plotFolder
    #
    # Silly trick to make the package checker happy... :-(
    lon <- lat <- Latitude <- Longitude <- label <- NULL
    #
    # Now compute the aggregation units, the values to be plotted, and make the map
    # Get the aggregated data for the aggregation units
    #
    aggLevel <- as.integer(map$aggregation)     		#; print(aggLevel)   
    aggUnitData <- map$aggUnitData				#; print(aggUnitData)
    #
    # For sample count markers, the colour may be based on a different admin division level from the aggregation
    #
    if (mapType %in% c("sampleCount","location")) {
        colourAdmDivLevel    <- param.getParam ("map.markerColourAggLevel", params)	#; print(colourAdmDivLevel) # should be 0 or 1
        colourAdmDivTitle    <- ADM_DIV_LABELS[colourAdmDivLevel+1]			#; print(colourAdmDivCol)
        colourAdmDivCol      <- GID_COLUMNS[colourAdmDivLevel+1]			#; print(colourAdmDivCol)
        colourAdmDivs        <- aggUnitData[,colourAdmDivCol]				#; print(colourAdmDivs)
        colourGids           <- unique(colourAdmDivs)					#; print(colourGids)
        colPalette           <- graphics.getColourPalette (userCtx)
        admDivPalette        <- rep_len(colPalette, length.out=length(colourGids))
        admDivTextPalette    <- graphics.makeTextPalette (admDivPalette)		
        names(admDivPalette) <- names(admDivTextPalette) <- colourGids			#; print(admDivPalette); print(admDivTextPalette)
        admDivPaletteLabels  <- map.getAdmDivNamesFromMeta (colourAdmDivLevel, sampleMeta)	#; print(admDivPaletteLabels)
    }
    #
    # Select the aggregation units to be plotted (all, in the case of location markers)
    # In this case, those that have an NA for the measure being plotted (should be only for drug resistance)
    #
    selAggUnitData <- aggUnitData				#; print(nrow(selAggUnitData))
    if (mapType != "location") {
        vals <- aggUnitData[,measureName]			#; print(vals)
        selAggUnitData <- aggUnitData[which(!is.na(vals)),]	#; print(nrow(selAggUnitData))
        if (nrow(selAggUnitData)==0) {
	    print(paste("Insufficient data for a map of", measureName, "for interval", interval$name))
	    return()
	}
    }
    #
    # Compute marker sizes. 
    #
    pointSizes <- markerMap.getAggUnitMarkerSizes (selAggUnitData, params)	#; print(pointSizes)
    #
    # Do the actual plot, starting with the background map
    #
    mapPlot <- baseMapInfo$baseMap
    #
    # If we need to show aggregation unit names, we need to compute the label positioning and plot
    #
    showMarkerNames <- param.getParam ("map.markerNames", params)
    if (showMarkerNames) {
        mapPlot <- map.addAggregationUnitNameLayer (mapPlot, selAggUnitData, baseMapInfo, params, markerSize=pointSizes)
    }
    #
    # Get the values to be displyed in the markers (except for the location markers)
    #
    if (mapType != "location") {
        mValues <- selAggUnitData[,measureName]
        valueLabels <- as.character(round(mValues, digits=2))
    }									#; print(valueLabels)
    #
    # Now add the markers, coloured according to the appropriate scale, depending on the type of map
    #
    if (mapType=="diversity") {
        markerColours <- param.getParam ("map.markerColours", params)
        # Two marker colours can be specified to create a gradient. If a single marker colour is 
        # specified, then create a gradient from white to that colour.
        if (length(markerColours) == 1) {
                    markerColours <- c("white", markerColours)
        }
        mMax <- max(mValues); scaleMax <- round(mMax,digits=1); scaleMax <- ifelse(scaleMax<mMax, scaleMax+0.1, scaleMax)
        mMin <- min(mValues); scaleMin <- round(mMin,digits=1); scaleMin <- ifelse(scaleMin>mMin, scaleMin-0.1, scaleMin)  #; print(paste(scaleMin, scaleMax))
        mapPlot <- mapPlot +
            ggplot2::geom_point(ggplot2::aes(x=Longitude, y=Latitude, fill=!!rlang::sym(measureName)),
                                data=selAggUnitData,
                                size=pointSizes, shape=21, stroke=params$markerBorderWidth) +
            ggplot2::scale_fill_gradientn(limits=c(scaleMin,scaleMax), colours=markerColours, values=c(0,1))
    } else if (mapType=="sampleCount") {			#; print(length(pointSizes)); print(params$markerBorderWidth)
        mapPlot <- mapPlot +
            ggplot2::geom_point(ggplot2::aes(x=Longitude, y=Latitude, fill=!!rlang::sym(colourAdmDivCol)),
	                        data=selAggUnitData,
                                size=pointSizes, shape=21, stroke=params$markerBorderWidth) +
            ggplot2::scale_fill_manual(values=admDivPalette, labels=admDivPaletteLabels, name=colourAdmDivTitle,
                                       guide=ggplot2::guide_legend(override.aes=list(
                                       size=params$legendFontSize, stroke=params$legendBorderWidth
                                       )))
    }  else if (mapType=="location") {
        mapPlot <- mapPlot +
            ggplot2::geom_point(ggplot2::aes(x=Longitude, y=Latitude, fill=!!rlang::sym(colourAdmDivCol)),
	                        data=selAggUnitData,
                                size=pointSizes, shape=21, stroke=params$markerBorderWidth) +
            ggplot2::scale_fill_manual(values=admDivPalette, labels=admDivPaletteLabels, name=colourAdmDivTitle,
                                       guide=ggplot2::guide_legend(override.aes=list(
                                           size=params$legendFontSize, stroke=params$legendBorderWidth
                                       )))
    } else if (mapType=="drug") {
        mapPlot <- mapPlot +
            ggplot2::geom_point(ggplot2::aes(x=Longitude, y=Latitude, fill=!!rlang::sym(measureName)), 
	                        data=selAggUnitData, 
	                        size=pointSizes, shape=21, stroke=params$markerBorderWidth) +
            ggplot2::scale_fill_gradientn(limits=c(0,1), colours=c("green3","orange2","red3","red3"), values=c(0, 0.2, 0.75, 1))
    } else if (mapType=="mutation") {
        markerColours <- param.getParam ("map.markerColours", params)
        #
        # Two marker colours can be specified to create a gradient. If a single marker colour is 
        # specified, then create a gradient from white to that colour.
        #
        if (length(markerColours) == 1) {
            markerColours <- c("white", markerColours)
        }
        scaleMin <- 0; scaleMax <- 1
        mapPlot <- mapPlot +
            ggplot2::geom_point(ggplot2::aes(x=Longitude, y=Latitude, fill=!!rlang::sym(measureName)), 
	                        data=selAggUnitData, size=pointSizes, shape=21, stroke=params$markerBorderWidth) +
            ggplot2::scale_fill_gradientn(limits=c(scaleMin,scaleMax), colours=markerColours, values=c(0,1))
    }
    #	    
    # Show the values in the markers
    #
    if (mapType=="sampleCount") {
        mapPlot <- mapPlot +
            ggplot2::geom_text(ggplot2::aes(x=Longitude, y=Latitude, colour=!!rlang::sym(colourAdmDivCol)),
                               data=selAggUnitData, 
                               label=valueLabels, hjust=0.5, vjust=0.5, size=params$map.markerValueFontSize, fontface="bold", 
                               show.legend=FALSE) +
            ggplot2::scale_colour_manual(values=admDivTextPalette)
    } else if (mapType != "location") {
        mapPlot <- mapPlot +
            ggplot2::geom_text(ggplot2::aes(x=Longitude, y=Latitude),
                               data=selAggUnitData,
                               label=valueLabels, hjust=0.5, vjust=0.5, size=params$map.markerValueFontSize, fontface="bold")
    }
    #
    # Return map plot for completion and saving to file
    #
    mapPlot
}
#
###############################################################################
# Graphical rendition of site markers
################################################################################
#
markerMap.getAggUnitMarkerSizes <- function(aggUnitData, params) {
    # Compute marker sizes. 
    # If only one size was given in the config, then the markers will be constant size/
    # If there are two sizes, then markers will be sized proportional to the number of samples, with the smaller
    # size representing 1 sample, and the larger size representing the numer of samples in the largest aggregation
    mSizeParam <- param.getParam ("map.markerSize", params)	#; print(mSizeParam) #; print(length(mSizeParam))
    if (length(mSizeParam) > 1) {
        minSize  <- mSizeParam[1]
        maxSize  <- mSizeParam[2]
        sizeRange  <- maxSize - minSize
        counts <- aggUnitData$SampleCount			#; print(counts)
        minCount <- 1
        maxCount <- max(counts)
        countRange <- maxCount - minCount
        markerSizes <- minSize + (sizeRange * ((counts - minCount) / countRange))
    } else {
        markerSizes <- mSizeParam
    }								#; print(markerSizes)
    markerSizes
}
#
###############################################################################
# Estimation of marker measures (drug resistance and diversity)
################################################################################
#
markerMap.resolveMeasureNames <- function(ctx, mapType, userMeasureNames) {
    measureNames <- userMeasureNames
    if (mapType=="diversity") {
        if ("ALL" %in% measureNames) {
            measureNames <- diversity.getDiversityMeasures()				#; print (measureNames)
        }
    } else if (mapType=="drug") {
        if ("ALL" %in% measureNames) {
            config <- context.getConfig(ctx)						#; print(names(config))
            measureNames <- setup.getFeatureNames(config$drugPredictionFeatures)	#; print (measureNames)
        }
    } else if (mapType=="mutation") {
        if ("ALL" %in% measureNames) {
            config <- context.getConfig(ctx)						#; print(names(config))
            measureNames <- setup.getFeatureNames(config$drugMutationFeatures)		#; print (measureNames)
        }
    }
    measureNames
}
#
markerMap.estimateMeasures <- function(ctx, datasetName, aggLevel, aggUnitData, sampleSetName, mapType, measureNames, params, dataFolder) {	#;print(measureNames)
    sampleMeta <- context.getMeta (ctx, datasetName)
    dataset <- ctx[[datasetName]]
    barcodeData  <- dataset$barcodes
    # Create aggregation index for each sample (the id of the aggregation unit where the sample originates)
    aggIndex <- map.getAggregationUnitIds (aggLevel, sampleMeta)
    
    # Get all aggregation units
    aggUnitGids <- rownames(aggUnitData)				#; print(aggUnitGids)
    
    # Get the data for all aggregation units
    measureData <- matrix(nrow=0, ncol=length(measureNames), dimnames=list(c(), measureNames))
    for (aIdx in 1:length(aggUnitGids)) {

        # Get the sample data to be aggregated for this unit
        aggUnitGid <- aggUnitGids[aIdx]					#; print(aggUnitGid)
        aggSamplesMeta <- sampleMeta[which(aggIndex == aggUnitGid),]	#; print(nrow(aggSamplesMeta))
        aggSamples <- rownames(aggSamplesMeta)
        
        # Get the admin division values from the first sample of this unit (assuming the values are the same for all)
        suffix <- ""
        if (mapType=="diversity") {
            cValues <- diversity.estimateDiversityMeasures (ctx, aggSamples, measureNames)
            suffix <- paste0("-", measureNames)
        } else if (mapType=="drug") {
            cValues <- meta.getResistancePrevalence (ctx, aggSamplesMeta, measureNames, params)
        } else if (mapType=="mutation") {
            cValues <- meta.getMutationPrevalence (ctx, aggSamplesMeta, measureNames, params)
        }                                                           	#; print(cValues)
        measureData <- rbind(measureData, cValues)                 	#; print(measureData)
    }
    measureData <- as.data.frame(measureData)				#;  print(measureData)
    aggUnitData <- cbind(aggUnitData, measureData)			#;  print(dim(aggUnitData))
    for (mIdx in 1:length(measureNames)) {
         measureName <- measureNames[mIdx]				#;  print(measureName)
         aggUnitData[,measureName] <- as.numeric(aggUnitData[,measureName])
    }

    # Write out the aggregation unit data to file
    aggDataFilename  <- paste0(dataFolder, "/AggregatedData-", sampleSetName, "-", aggLevel, suffix, ".tab")
    utils::write.table(aggUnitData, file=aggDataFilename, sep="\t", quote=FALSE, row.names=FALSE)

    aggUnitData
}
