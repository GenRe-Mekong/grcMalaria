###############################################################################
# Cluster Analysis 
###############################################################################
#
# Main entry point- load stored clusters, or identify them from scratch if they are not stored
#
cluster.findClusters <- function (userCtx, sampleSetName, params) {
    sampleSet <- context.getSampleSet (userCtx, sampleSetName)
    if (is.null(sampleSet)) {
        stop(paste("Sample set not initialized:", sampleSetName))
    }
    ctx       <- sampleSet$ctx
    config    <- context.getConfig(ctx) 

    clusterSetName    <- param.getParam ("cluster.clusterSet.name", params)	#; print(clusterSetName)
    minIdentityLevels <- param.getParam ("cluster.identity.min", params) 	#; print(minIdentityLevels)
    useImputation     <- param.getParam ("cluster.impute", params)		#; print(useImputation)
    method            <- param.getParam ("cluster.method", params)		#; print(method)
    print(paste("Clustering Method:",method))
    
    dataRootFolder  <- getOutFolder(config, sampleSetName, c("cluster", "data", clusterSetName), create=FALSE)
    if (file.exists(dataRootFolder)) {
        unlink(dataRootFolder, recursive=TRUE)
    }

    clusterSetInfos <- list()
    for (idIdx in 1:length(minIdentityLevels)) {
        minIdentity <- minIdentityLevels[idIdx]
        minIdentityLabel <- getMinIdentityLabel (minIdentity)			#; print(minIdentityLabel)

        # Determine what clustering approach must be used.
        # TODO - At the moment, only graph-based clustering methods are implemented, but this may be extended later.
        clustersData <- NULL
        if (method %in% cluster.graphMethods) {
            clustersData <- cluster.findClustersFromGraph (ctx, sampleSetName, clusterSetName, method, minIdentity, useImputation, params)
        } else {
            stop (paste("Invalid method specified in parameter 'cluster.method':", method))
        }
    
        # Sort clusters by descending size, and label them in order, appending the sequence number to the cluster set name
        clustersData <- clustersData[rev(order(clustersData$Count)),]		#; print(head(clustersData))
        clusterSetName <- param.getParam ("cluster.clusterSet.name", params)	#; print(clusterSetName)
        clustersData$ClusterId <- paste(clusterSetName, 
                                        formatC(seq(1,nrow(clustersData)), width=3, format="d", flag="0"), sep="-")
        rownames(clustersData) <- clustersData$ClusterId			#; print(head(clustersData))
        clustersData <- as.data.frame(clustersData[,c("ClusterId","Count","SampleList")])

        # Write out result files
        dataFileSuffix <- paste("", sampleSetName, clusterSetName, minIdentityLabel, sep="-")
        dataFolder  <- getOutFolder(config, sampleSetName, 
                                    c("cluster", "data", clusterSetName, minIdentityLabel))

        # Write out cluster definitions
        clustersDataFile <- paste0(dataFolder, "/clusters", dataFileSuffix, ".tab")
        utils::write.table(clustersData, file=clustersDataFile, quote=TRUE, sep="\t", row.names=FALSE, col.names=TRUE)

        # Write out sample/cluster association
        memberData <- cluster.getMemberData (clustersData)
        memberDataFile <- paste0(dataFolder, "/clusterMembers", dataFileSuffix)
        writeSampleData(memberData, memberDataFile)

        # Write out cluster stats
        statsData <- cluster.getClusterStats (ctx, clustersData, memberData)
        statsDataFile <- paste0(dataFolder, "/clusterStats", dataFileSuffix)
        writeSampleData(statsData, statsDataFile)

        # Keep the cluster info for storing in the context
        clusterSetInfo <- list(clusterSetName=clusterSetName,
                               sampleSetName=sampleSetName,
                               minIdentity=minIdentity,
                               clusters=clustersData,	# dataframe
                               members=memberData, 	# dataframe
                               stats=statsData)		# dataframe
        clusterSetInfos[[minIdentityLabel]] <- clusterSetInfo
    }

    # Reference the cluster data from the context
    sampleSet$clusters[[clusterSetName]] <- clusterSetInfos	#; print(names(sampleSet$clusters))
}
#
# Retrieve the cluster data from the context
#
cluster.getClustersSetFromContext <- function(userCtx, sampleSetName, clusterSetName) {
    sampleSet <- context.getSampleSet (userCtx, sampleSetName)	#; print(names(userCtx)); print(names(userCtx$sampleSets))
    if (is.null(sampleSet)) {
        stop(paste("Sample set not found:", sampleSetName))        
    }								#; print(names(sampleSet));    
    clusterSetInfos <- sampleSet$clusters[[clusterSetName]]
    if (is.null(clusterSetInfos)) {
       stop(paste("Invalid cluster set specified:", clusterSetName))        
    }
    clusterSetInfos
}
#
#
#
cluster.getMemberData <- function(clustersData) {
    cIds <- c()
    sIds <- c()
    for (clIdx in 1 : nrow(clustersData)) {
         # Get members of the cluster and label them
         clusterName <- clustersData$ClusterId[clIdx]
         clSampleNames <- unlist(strsplit(clustersData$SampleList[clIdx], split=","))
         sIds <- c(sIds, clSampleNames)
         cIds <- c(cIds, rep(clusterName, length(clSampleNames)))
    }
    #print(sIds)
    #print(cIds)
    clusterMembers <- data.frame(Sample=sIds, Cluster=cIds, stringsAsFactors=FALSE)
    rownames(clusterMembers) <- sIds
    clusterMembers
}

#
#
#
cluster.getClusterPalette <- function(ctx, clusterIds) {
    # Get the name of all the clusters and sort them 
    clusterIds <- sort(clusterIds)
    # Put "Other" at the end if it is there
    if ("Other" %in% clusterIds) {
        clusterIds <- c(clusterIds[which(clusterIds != "Other")], "Other")
    }									#; print(clusterIds)
    # Construct a palette using the default palette, recycling it if there are too many clusters
    # and adding white as the last colour for class "Other" (samples not belionging to a cluster)
    colPalette <- graphics.getColourPalette (ctx)
    clusterPalette <- rep_len(colPalette, length.out=(length(clusterIds)-1))
    clusterPalette <- c(clusterPalette,"white")
    names(clusterPalette) <- clusterIds					#; print(clusterPalette)
    clusterPalette
}
#
# #######################################################################################
#
# Descriptive data about the clusters (e.g. prevalence of mutations, etc.)
#
cluster.getClusterStats <- function(ctx, clustersData, clusterMembers) {
    config <- context.getConfig(ctx)

    clNames <- rownames(clustersData)			#; print(head(clustersData)); print(clNames)
    
    # Create a table of stats data
    statsData  <- NULL
    statsNames <- NULL
    for (clIdx in 1 : length(clNames)) {
        clName <- clNames[clIdx]			#; print(clName)
        
        # Get the metadata for this cluster
        clSampleNames <- clusterMembers$Sample[which(clusterMembers$Cluster==clName)]
        clSampleMeta <- context.getSampleMeta (ctx, clSampleNames)
        
        # Get the sample count first
        statsNames <- "Count"
        statValues <- length(clSampleNames)

        # Get the drug resistance prevalences for this cluster
        drugNames  <- setup.getFeatureNames(config$cluster.stats.drugPredictionFeatures)
        if (!is.null(drugNames)) { 
            clPrevalence <- meta.getResistancePrevalence (ctx, clSampleMeta, drugNames)
            clPrevalence <- format(as.numeric(clPrevalence), digits=2, nsmall=2)
            statsNames  <- c(statsNames, drugNames)
            statValues <- c(statValues, clPrevalence)
        }

        # Get the counts for this cluster
        countableFeatureNames <- setup.getFeatureNames(config$cluster.stats.countableFeatures)
        if (!is.null(countableFeatureNames)) {
            clCounts <- meta.getValueCounts (ctx, clSampleMeta, countableFeatureNames)
            statsNames <- c(statsNames, countableFeatureNames)
            statValues <- c(statValues, clCounts)
        }

        # Get the mutation prevalences for this cluster
        mutationNames  <- setup.getFeatureNames(config$cluster.stats.drugMutationFeatures)
        if (!is.null(mutationNames)) { 
            clPrevalence <- meta.getMutationPrevalence (ctx, clSampleMeta, mutationNames, params=NULL)
            clPrevalence <- format(as.numeric(clPrevalence), digits=2, nsmall=2)
            statsNames <- c(statsNames, mutationNames)
            statValues <- c(statValues, clPrevalence)
        }

        # Stick the joined row data to the Stats table
        names(statValues) <- statsNames
        statsData <- rbind(statsData, statValues)
    }
    statsData <- data.frame(statsData)
    colnames(statsData) <- statsNames
    rownames(statsData) <- clNames
    statsData
}
#
cluster.getClusterStatsText <- function(clusterStats, clustersName) {
    statsNames  <- colnames(clusterStats)
    statsValues <- clusterStats[clustersName,]
    clusterInfoTextLines <- paste(statsNames, statsValues, sep=": ")
    clusterInfoText <- paste(clusterInfoTextLines, collapse = "\n")		#; print(clusterInfoText)
    clusterInfoText
}
#
###############################################################################
# Graph-based Clustering
###############################################################################
cluster.graphCommunityMethods <- c("louvain", "leiden")
cluster.graphMethods          <- c("allNeighbours", cluster.graphCommunityMethods)

cluster.findClustersFromGraph <- function (ctx, sampleSetName, clusterSetName, method, minIdentity, useImputation, params) {
    config <- context.getConfig(ctx)
    
    # Get a table of pairwise distance/identity values for all pairs of samples that meet the threshold
    distData  <- context.getDistanceMatrix (ctx, sampleSetName, useImputation)	#; print(head(distData))#; print(nrow(distData))

    edgeData <- clusterGraph.getPairwiseIdentityData (distData, minIdentity, params)	#; print(head(edgeData))
    edgeData$weight <- (edgeData$Identity * edgeData$Identity)			#; print(nrow(edgeData))
    
    nodeNames <- unique(c(as.character(edgeData$Sample1),as.character(edgeData$Sample2)))
    nodeCount <- length(nodeNames)						#; print (paste(length(nodeNames),length(unique(edgeData$Sample1)),length(unique(edgeData$Sample2))))
    nodeData  <- data.frame(NodeName=nodeNames, Count=rep(1,nodeCount), NodeType=rep("sample",nodeCount))
    gr <- igraph::graph_from_data_frame(edgeData, directed=FALSE, vertices=nodeData)	#; print(paste("processClusters",minIdentity))

    # Perform clustering from the graph, identifying all clusters of sufficient size
    if (method == "allNeighbours") {
        clustersList <- cluster.findAllNeighbourClusters (gr, params)
    } else if (method %in% cluster.graphCommunityMethods) {
        clustersList <- cluster.findGraphCommunities (gr, method, params)
    }
    
    # Check if we actually found any clusters
    clCount <- length(clustersList)						#; print(clCount)
    if (clCount == 0) {
        print(paste("No clusters of desired minimum side were found at identity threshold", minIdentity))
        return (NULL)
    }
    
    # Turn the list of clusters into a data frame    
    sampleCounts <- integer(clCount)
    sampleLists <- character(clCount)
    for (clIdx in 1 : clCount) {
        clSampleNames <- clustersList[[clIdx]]
        sampleCounts[clIdx] <- length(clSampleNames)
        sampleLists[clIdx] <- paste(clSampleNames, collapse=",")
    }
    clustersData <- data.frame(Count=sampleCounts, SampleList=sampleLists, 
                               stringsAsFactors=FALSE)				#; print(head(clustersData))
    clustersData
}
#
# #######################################################################################
#
# Clustering Method: identification by traversing All Neighbours.
#
# To find each cluster, we pick a node, and find all neighbouts recursively.
#
cluster.findAllNeighbourClusters <- function (gr, params) {
    # Make a working copy of the graph and count the nodes
    wg <- igraph::induced_subgraph(gr, igraph::V(gr))	#; plot(wg)
    nodeCnt <- length(igraph::V(wg))    		#; print(paste("nodeCnt",nodeCnt))

    # Identify all clusters with >= minCount samples
    minCount <- param.getParam ("cluster.minSize", params)
    clList <- list()
    while (nodeCnt > 0) {				#; print(nodeCnt)
        #
        # Get the fist node in the graph, and build a cluster by traversing to all neighbours
        #
        firstNode <- igraph::V(wg)[1]
        clNodeIdxs <- cluster.findAllNeighbourConnectedNodes (wg, firstNode, as.integer(firstNode))
        						#; print(paste("new cluster - nodes:",paste(clNodeIdxs,collapse=",")))
        clSampleCount <- length(clNodeIdxs)		#; print(clSampleCount)
        if (clSampleCount >= minCount) {
            clSampleNames <- names(igraph::V(wg))[clNodeIdxs]
            clIdx <- length(clList)+1
            clList[[clIdx]] <- clSampleNames
        }
        # Remove the samples from the cluster we found
        keepSamples <- seq(1, nodeCnt)[-clNodeIdxs]
        wg <- igraph::induced_subgraph(wg, keepSamples, impl="copy_and_delete")
        nodeCnt <- length(igraph::V(wg))
    }
    clList
}

cluster.findAllNeighbourConnectedNodes <- function (gr, curr, nodes) {
    ns <- igraph::neighbors(gr, curr)
    if (length(ns) > 0) {
        foundNew <- FALSE
        for (nsIdx in 1 : length(ns)) {
            n <- ns[nsIdx]
            nl <- as.integer(n)
            if (!(nl %in% nodes)) {
                foundNew <- TRUE
                nodes <- c (nodes, nl)
                nodes <- cluster.findAllNeighbourConnectedNodes (gr, n, nodes)
            }
        }
    }
    nodes
}
#
# #######################################################################################
#
# Clustering Method: Community Analysis.
#
cluster.findGraphCommunities <- function (gr, method, params) {

    if (method == "louvain") {
        resolution <- param.getParam ("cluster.resolution", params)
        partition <- igraph::cluster_louvain(gr, weights=igraph::E(gr)$weight, resolution_parameter=resolution) 

    } else if (method == "leiden") {
        objFunction <- param.getParam ("cluster.objFunction", params)
        resolution  <- param.getParam ("cluster.resolution", params)
        beta        <- param.getParam ("cluster.beta", params)
        partition <- igraph::cluster_leiden(gr, objective_function=objFunction, weights=igraph::E(gr)$weight, resolution_parameter=resolution, beta=beta) 

    } else {
        stop (paste("Invalid method specified in parameter 'cluster.method':", method))
    }
    nodeComms <- partition$membership
    commIds <- sort(as.integer(unique(nodeComms)))
    sampleNames <- names(igraph::V(gr))
    #names(nodeComms) <- sampleNames
    
    minCount <- param.getParam ("cluster.minSize", params)
    clList <- list()
    for (clIdx in 1:length(commIds)) {
        commId <- commIds[clIdx]
        clNodeIdxs <- which(nodeComms == commId)
        clSampleCount <- length(clNodeIdxs)		#; print(clSampleCount)
        if (clSampleCount < minCount) {
            next
        }
        clSampleNames <- sampleNames[clNodeIdxs]
        clIdx <- length(clList)+1
        clList[[clIdx]] <- clSampleNames
    }
    clList
}
