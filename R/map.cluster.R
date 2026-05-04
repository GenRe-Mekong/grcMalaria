###############################################################################
# Map Cluster Analysis
################################################################################
#
CLUSTER_MAP_PREFIX <- "clusterPrevalence-ge"
#
################################################################################
#
clusterMap.executeMap <- function (map) {

    mapMaster   <- map$master
    mapType     <- mapMaster$type
    measureName <- map$measureName						#; print(measureName)
    interval    <- map$interval
    #
    datasetName <- map$datasetName
    sampleSet   <- mapMaster$sampleSet
    userCtx     <- mapMaster$userCtx
    params      <- mapMaster$params
    config      <- context.getConfig(userCtx)
    #
    sampleMeta <- context.getMeta (map$mapCtx, datasetName) 
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
    lon <- lat <- label <- x <- y <- alpha <- NULL
    #
    # Parse the measure to get the cluster name and the minIdentity level
    #
    mParts <- clusterMap.parseMeasureName (measureName)				#; print(measureName)
    minIdentity <- mParts$minIdentity						#; print(minIdentity)
    clusterId <- mParts$clusterId						#; print(clusterId)		
    #
    # Now compute the aggregation units, the values to be plotted, and make the map
    # Get the aggregated data for the aggregation units
    #
    aggLevel <- as.integer(map$aggregation)     				#; print(aggLevel)   
    aggUnitData <- map$aggUnitData						#; print(aggUnitData)
    #    
    # Get the cluster data
    #
    clusterSetLabel  <- getMinIdentityLabel (minIdentity)			#; print(clusterSetLabel)
    clusterSet       <- mapMaster$clusterSets[[clusterSetLabel]]		#; print(str(clusterSet))
    #
    clusterShareData <- clusterMap.getAggUnitClusterCountData (clusterSet, sampleMeta, aggLevel, aggUnitData, params)
											#; print(clusterShareData)
    clusterPalette <- mapMaster$clusterSetPalettes[[clusterSetLabel]]
    clusterInfoText <- cluster.getClusterStatsText(clusterSet$stats, clusterId)		#; print(clusterInfoText)
    #
    #
    #
    clusterShareData <- clusterShareData[which(clusterShareData$Cluster == clusterId),]	#; print(head(clusterShareData))
    selAggUnitIds <- as.character(clusterShareData$UnitId)
    selAggUnitData <- aggUnitData[selAggUnitIds,]					#; print(nrow(selAggUnitData))
    #
    # Do the actual plot, starting with the background map
    #
    mapPlot <- baseMapInfo$baseMap
    #
    #
    #
    mapPlot <- clusterMap.addConnections (mapPlot, baseMapInfo, clusterId, clusterShareData, clusterPalette, params)
    #        
    # Dirty trick so we can show an annotation with cluster info above the legends.
    # We "plot" a couple of points, and create a ficticious alpha legend containing the cluster info text.
    #
    bb <- baseMapInfo$anBB
    dummydf <- data.frame(x=c(bb$xMin,bb$xMin), y=c(bb$yMax,bb$yMax), alpha=c(0.1, 0.11))
    mapPlot <- mapPlot +
        ggplot2::geom_point(ggplot2::aes(x=x, y=y, alpha=as.numeric(alpha)), size=0.01, 
                            data=dummydf) +
        ggplot2::scale_alpha_continuous(clusterId, breaks=c(0.1, 0.11), labels=c(clusterInfoText,""), 
                                        guide=ggplot2::guide_legend(order=1,keywidth=0, keyheight=0.01, nrow=1,
                                        override.aes=list(shape=NA,fill=NA,size=0.01)))
    #
    plotTitle <- paste(sampleSet$name, "-", clusterId);
    mapPlot <- mapPlot + ggplot2::labs(title=plotTitle, subtitle="")
    #
    # Return map plot for completion and saving to file
    #
    mapPlot
}
#
################################################################################
#
clusterMap.resolveMeasureNames <- function(clusterSets) {
    expanded <- c()
    minIdentityLabels <- names(clusterSets)
    setCount <- length(minIdentityLabels)
    for (idIdx in 1:setCount) {
        minIdentityLabel <- minIdentityLabels[idIdx]
        minIdentity <- getMinIdentityFromLabel (minIdentityLabel)
        clusterSet <- clusterSets[[minIdentityLabel]]
        clusterData <- clusterSet$clusters
        clusterIds <- as.character(clusterData$ClusterId)
        measureLabels <- paste(getMinIdentityLabel(minIdentity, prefix=CLUSTER_MAP_PREFIX), 
                               clusterIds, sep="/")				#; print(measureLabels)
        expanded <- c(expanded, measureLabels)
    }										#; print(expanded)
    expanded
}
#
################################################################################
#
clusterMap.parseMeasureName <- function(measureName) {
    prefix <- CLUSTER_MAP_PREFIX
    if (!startsWith(measureName, prefix)) {
        return (NULL)				#; print("Incorrect prefix")
    }
    suffix <- substring(measureName,nchar(prefix)+1)
    sParts <- unlist(strsplit(suffix, "/"))
    list(minIdentity=as.numeric(sParts[1]), clusterId=sParts[2])
}
#
################################################################################
# 
# For each cluster set, construct a palette using the default palette, recycling colours if there are too 
# many clusters, adding white as the last colour for class "Other" (samples not belionging to a cluster)
#
clusterMap.getClusterSetsPalettes <- function (ctx, clusterSets) {
    palettes <- list()
    minIdentityLabels <- names(clusterSets)
    setCount <- length(minIdentityLabels)
    for (idIdx in 1:setCount) {
        minIdentityLabel <- minIdentityLabels[idIdx]
        clusterSet <- clusterSets[[minIdentityLabel]]
        clusterData <- clusterSet$clusters
        clusterIds <- as.character(clusterData$ClusterId)
	palettes[[minIdentityLabel]] <- cluster.getClusterPalette (ctx, clusterIds)
    }										#; print(expanded)
    palettes
}
#
################################################################################
#
#
#
clusterMap.getAggUnitClusterCountData <- function(clusterSet, sampleMeta, aggLevel, aggUnitData, params) {
    #
    # clusterData is a dataframe with three columns: ClusterId, Count (of samples), SampleList (comma-separated lists of sample ids)
    #
    clusterData <- clusterSet$clusters							#; print(head(clusterData))
    #
    # memberData is a dataframe with three columns: Sample, ClusterId. It is in essence a derivative of the clusterData
    #
    memberData  <- clusterSet$members							#; print(nrow(memberData)); print(head(memberData))
    #
    clusterData <- clusterData[order(clusterData$Cluster),]				#; print(nrow(clusterData))
    clusterIds <- as.character(clusterData$Cluster)					#; print(clusterIds)
    #
    # It is *possible* that some members of the clusters are *not* in the metadata (e.g. if imputation is used in cluster selection or not)
    #
    sampleNames <- rownames(sampleMeta)
    memberData <- memberData[which(memberData$Sample %in% sampleNames),]		#; print(length(sampleNames))
    memberIds      <- memberData$Sample							#; print(length(memberIds))
    memberClusters <- memberData$Cluster						#; print(memberClusters)
    #
    # Assign cluster Ids to the sample metadata
    #
    sampleCluster <- rep("-", length(sampleNames))
    names(sampleCluster) <- sampleNames
    sampleCluster[memberIds] <- memberClusters
    sampleMeta$Cluster <- sampleCluster
    #
    # Get all aggregation unit ids
    #
    aggUnitGids <- rownames(aggUnitData)						#; print(aggUnitGids); print(nrow(sampleMeta))
    #
    # Create aggregation index for each sample (the id of the aggregation unit where the sample originates)
    #
    aggIndex <- map.getAggregationUnitIds (aggLevel, sampleMeta)
    #
    # Get the data for all aggregation units
    #
    countData <- NULL
    for (aIdx in 1:length(aggUnitGids)) {
        # Get the sample data to be aggregated for this unit
        aggUnitGid <- aggUnitGids[aIdx]
        unitMeta <- sampleMeta[which(aggIndex == aggUnitGid),]				#; print(nrow(unitMeta))
        clUnitMeta <- unitMeta[which(unitMeta$Cluster != "-"),]
        unitGroupCounts <- table(clUnitMeta$Cluster)
        unitGroupCounts <- unitGroupCounts[order(-unitGroupCounts)]
        unitGroupNum <- length(unitGroupCounts)+1

        unitData <- aggUnitData[aIdx,]
        df <- data.frame(matrix(unitData,ncol=ncol(aggUnitData),nrow=unitGroupNum,byrow=TRUE))
        colnames(df) <- colnames(aggUnitData)
        df$Cluster <- c(names(unitGroupCounts), "Other")
        df$ClusterCount <- c(unitGroupCounts, nrow(unitMeta)-nrow(clUnitMeta))
        df$ClusterProp <- df$ClusterCount / unitData$SampleCount
        countData <- rbind(countData, df)
    }
    countData$Cluster <- factor(countData$Cluster)
    countData$Latitude <- as.numeric(countData$Latitude)
    countData$Longitude <- as.numeric(countData$Longitude)
    countData$SampleCount <- as.numeric(countData$SampleCount)				#; print(countData)
    countData						
}
#
###############################################################################
# Maps of connections between Haplotype Sharing sites plotting 
################################################################################
#
clusterMap.addConnections <- function (mapPlot, baseMapInfo, clusterId, clusterShareData, clusterPalette, params) { #; print(clusterId) ; print(clusterShareData)
    #
    # Work out the values to display for every marker (the fractional prevalence)
    #
    mValues <- as.numeric(clusterShareData$ClusterProp)			#; print(str(clusterShareData))
    valueLabels <- round(mValues, digits=2)
    
    # Work out size, colour and colour gradient for the markers
    pointSizes <- params$map.markerSize					#; print(clusterPalette)
    clusterColour <- clusterPalette[clusterId]				#; print(clusterColour)
    mMax <- max(mValues); scaleMax <- round(mMax,digits=1); scaleMax <- ifelse(scaleMax<mMax, scaleMax+0.1, scaleMax)
    mMin <- min(mValues); scaleMin <- round(mMin,digits=1); scaleMin <- ifelse(scaleMin>mMin, scaleMin-0.1, scaleMin)  #; print(paste(scaleMin, scaleMax))

    # For ech pair, work out a connectedness measure, in this case the fraction of sample pairs in which both samples carry this cluster     
    aggUnitCount <- nrow(clusterShareData)
    if (aggUnitCount > 1) {
        connectData <- NULL
        for (a1Idx in 1:(aggUnitCount-1)) {
            for (a2Idx in (a1Idx+1):aggUnitCount) {
                # Get the sample data
                lat1 <- clusterShareData$Latitude[a1Idx]
                lon1 <- clusterShareData$Longitude[a1Idx]
                sampleCount1 <- clusterShareData$SampleCount[a1Idx]
                clusterCount1 <- clusterShareData$ClusterCount[a1Idx]
                #
                lat2 <- clusterShareData$Latitude[a2Idx]
                lon2 <- clusterShareData$Longitude[a2Idx]
                sampleCount2 <- clusterShareData$SampleCount[a2Idx]
                clusterCount2 <- clusterShareData$ClusterCount[a2Idx]
                #
                pairCount <- sampleCount1 * sampleCount2
                matchCount <- clusterCount1 * clusterCount2
                prop <- matchCount / pairCount
            
                cValues <- c(lat1, lon1, lat2, lon2, prop)
                connectData <- rbind(connectData, cValues)
            }
        }
        aggUnitPairData <- data.frame(connectData)						#; print (dim(aggUnitPairData))
        nc <- ncol(aggUnitPairData)
        aggUnitPairData[,1:nc] <- sapply(aggUnitPairData[,1:nc], as.numeric)			#; print(ncol(aggUnitPairData))
        colnames(aggUnitPairData) <- c("Lat1", "Lon1", "Lat2", "Lon2", "Weight")		#; print(aggUnitPairData)
        
        hc <- grDevices::rgb2hsv(grDevices::col2rgb(clusterColour))
        weakHaploColour <- grDevices::hsv(h=hc[1,1],s=hc[2,1]*0.2,v=hc[3,1])
        
        # Silly trick to make the package checker happy... :-(
        Lon1 <- Lat1 <- Lon2 <- Lat2 <- Weight <- NULL
	
        # Now plot the connections
        #print(aggUnitPairData); #print(nrow(aggUnitPairData)); print(colnames(aggUnitPairData))
        mapPlot <- mapPlot +
                   ggplot2::geom_curve(ggplot2::aes(x=Lon1, y=Lat1, xend=Lon2, yend=Lat2, linewidth=Weight, colour=Weight),	# draw edges as arcs
                               data=aggUnitPairData, 
                               curvature=0.2, alpha=0.75) +
                   ggplot2::scale_linewidth_continuous(guide="none", range=c(params$connectCurveWidthMin, params$connectCurveWidthMax)) + 
                   ggplot2::scale_colour_gradientn(name="Sharing", colours=c(weakHaploColour,clusterColour))
    }
    #
    # Silly trick to make the package checker happy... :-(
    Longitude <- Latitude <- ClusterProp <- NULL
    #
    # If we need to show aggregation unit names, we need to compute the label positioning and plot
    #
    showMarkerNames <- param.getParam ("map.markerNames", params)
    if (showMarkerNames) {
        mapPlot <- map.addAggregationUnitNameLayer (mapPlot, clusterShareData, baseMapInfo, params, markerSize=pointSizes)
    }
    #
    # Plot the Aggregation Unit markers
    mapPlot <- mapPlot +
	       ggplot2::geom_point(ggplot2::aes(x=Longitude, y=Latitude, fill=ClusterProp), 
	                           data=clusterShareData, 
	                           size=pointSizes, shape=21, stroke=params$markerBorderWidth) +
	       ggplot2::scale_fill_gradientn(name="Prevalence", limits=c(scaleMin,scaleMax), colours=c("white",clusterColour), values=c(0,1)) +
               ggplot2::geom_text(ggplot2::aes(x=Longitude, y=Latitude), 
                                  data=clusterShareData, 
                                  label=valueLabels, size=params$map.markerValueFontSize, fontface="bold", hjust=0.5, vjust=0.5)
    mapPlot
}


