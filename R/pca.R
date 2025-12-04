phylo.pcCount <- 10
#
###############################################################################
# PCA Analysis.
# Computes and plots PCA, using a pairwise distance matrix as an input.
# Multiple analyses, based on different subsets of samples, can be performed in a single run.
#
# Supported methods:
#     pca/PCoA		Analyzes a distance matrix
################################################################################
#
pca.execute <- function(userCtx, sampleSetName, method, params) {
    sampleSet <- context.getSampleSet (userCtx, sampleSetName)
    ctx <- sampleSet$ctx
    #
    # Get the output folders
    #
    config <- context.getConfig(ctx)
    dataFolder      <- getOutFolder(config, sampleSetName, c(method, "data"))
    plotsRootFolder <- getOutFolder(config, sampleSetName, c(method, "plots"))
    #
    # Get the metadata, distance and genotypes.
    # For MCA we have to use imputed data to remove hets and missing
    #
    if (method == "MCA") {
        datasetName <- "imputed"
    } else {
        useImputed <- param.getParam ("impute", params)				#; print(useImputed)
        datasetName <- ifelse (useImputed, "imputed", "filtered")
    }
    sampleMeta <- context.getMeta (ctx, datasetName) 
    #
    # Compute the Principal components
    #
    pcNames <- paste("PC", seq(1:phylo.pcCount), sep="")
    sampleNames <- rownames(sampleMeta)
    #
    pcScores <- NULL
    varExplained <- NULL
    if (method == "PCoA") {
        distData  <- context.getDistanceMatrix (ctx, sampleSetName, useImputed)
        pcaResults <- stats::cmdscale(as.matrix(distData), eig=TRUE, k=phylo.pcCount)
        pcScores <- pcaResults$points
        varExplained <- pcaResults$eig / sum(abs(pcaResults$eig))

    } else if (method == "MCA") {
        catData <- ctx$rootCtx$imputed$barcodeGenoTable
        catData <- catData[sampleNames,]
        catData[] <- lapply(catData, as.factor)
        mcaResult <- FactoMineR::MCA(catData, ncp=10, graph=FALSE)
        pcScores <- mcaResult$ind$coord
        varExplained <- mcaResult$eig[,"percentage of variance"]
        
        varCoords <- data.frame(Allele=rownames(mcaResult$var$coord), mcaResult$var$coord)
        varCoordsFilename  <- paste(dataFolder, paste0(pca.getDataFileName("-variantCoords", sampleSetName, method),".tab"), sep="/")
        utils::write.table(varCoords, file=varCoordsFilename, sep="\t", quote=FALSE, row.names=FALSE)
    }
    #
    # Write out the data to file
    #
    rownames(pcScores) <- sampleNames
    colnames(pcScores) <- pcNames
    pcaScoresFilename  <- paste(dataFolder, pca.getDataFileName("", sampleSetName, method), sep="/")
    writeSampleData(pcScores, pcaScoresFilename)
    #
    pcaVarData <- data.frame(pcNames, varExplained[1:length(pcNames)])
    colnames(pcaVarData) <- c("PC","VarianceExplained")
    pcaVarFilename  <- paste(dataFolder, paste0(pca.getDataFileName("-varExplained", sampleSetName, method),".tab"), sep="/")
    utils::write.table(pcaVarData, file=pcaVarFilename, sep="\t", quote=FALSE, row.names=FALSE)
    #
    # Attach the PCA results to the metadata
    #
    sampleMeta <- cbind(sampleMeta, pcScores)
    #
    # Get the plot definitions and execute them
    #
    plotList <- param.getParam ("plot.plotList", params)		#; print(plotList)
    for (plotIdx in 1:length(plotList)) {
        plotDef <- plotList[[plotIdx]]
        plotDefName <- plotDef$name
        plotName <- paste(sampleSetName, plotDefName, method, sep="-")
        print (paste("PCA Plot: ",plotName))
        #
        # Set up the graphical attributes for rendering
        #
        ga <- graphics.getGraphicalAttributes(userCtx, sampleMeta, plotDef$attributes)
        gaData <- ga$sampleAttrData
        legendData <- ga$legendData
        #
        # Merge graphics attributes and metadata
        #
        plotMetadata <- cbind(sampleMeta, gaData)
        #
        # Write out the plot data (facilitate debug)
        #
        plotMetadataFilename  <- paste0(dataFolder, "/", "plotData-", plotName, ".tab")
        writeSampleData(plotMetadata, plotMetadataFilename)
        #
        # Order the samples so they are plotted in the correct stack order
        #
        excludedIdx <- which(plotMetadata$plot__order <= 0)        	#; print(plotMetadata$plot__order)
        if (length(excludedIdx) > 0) {
            print(paste("Excluded samples:", length(excludedIdx)))
            plotMetadata <- plotMetadata[-excludedIdx,]
        }
        plotMetadata <- plotMetadata[order(-plotMetadata$plot__order),]

        # Do the plots
        plotsFolder  <- getSubFolder (plotsRootFolder, plotDefName)
        pca.plotPrincipalComponents (plotName, plotMetadata, legendData, "PC1", "PC2", "PC3", plotsFolder, params);
        pca.plotPrincipalComponents (plotName, plotMetadata, legendData, "PC1", "PC4", "PC5", plotsFolder, params);
        pca.plotPrincipalComponents (plotName, plotMetadata, legendData, "PC1", "PC6", "PC7", plotsFolder, params);
        pca.plotPrincipalComponents (plotName, plotMetadata, legendData, "PC1", "PC8", "PC9", plotsFolder, params);
    }
}

pca.getDataFileName <- function(suffix, sampleSetName, method) {
    fn <- paste("pca", suffix, "-", sampleSetName, "-", method, sep="")
    fn
}

###############################################################################
# Principal Component plotting
###############################################################################
pca.plotPrincipalComponents <- function(plotName, sampleMeta, legendData, 
                                    pcAIdx, pcBIdx, pcCIdx=NULL, plotsFolder, params) {
    plotThree <- !is.null(pcCIdx)
    sampleCount <- nrow(sampleMeta)
  
    #Get the graphical parameters for each sample
    sampleColours <- sampleMeta$plot__colour		#; print(sampleColours)
    samplePch     <- sampleMeta$plot__pch		#; print(samplePch)
    sampleSize    <- sampleMeta$plot__size		#; print(sampleSize)
    sampleLwd     <- sampleMeta$plot__lwd		#; print(sampleLwd)
    sampleLcolours<- sampleMeta$plot__lcolour		#; print(sampleLcolours)
  
    # Plot
    plotFilename  <- paste("pca-",plotName,"-",pcAIdx,"_",pcBIdx, sep="")
    if (plotThree) {
        plotFilename  <- paste(plotFilename, pcCIdx, sep="_")
    }
    print(paste("Plotting ", plotFilename))
    graphicFilenameRoot  <- paste(plotsFolder, plotFilename, sep="/")
    initializeGraphics (getGraphicsFilename (graphicFilenameRoot), params)
    if (plotThree) {
        graphics::par(mfrow=c(2,2))
    } else {
#        graphics::par(mfrow=c(1,2))
    }
    graphics::par(mar=c(4,4,1,1) + 0.1)
    x <- sampleMeta[,pcAIdx]
    y <- sampleMeta[,pcBIdx]
    plot(x,y, xlab=pcAIdx, ylab=pcBIdx, cex.lab=0.8, cex.axis=0.8, bg=sampleColours, col=sampleLcolours, pch=samplePch, lwd=sampleLwd, cex=sampleSize)
    if (plotThree) {
        x <- sampleMeta[,pcCIdx]
        y <- sampleMeta[,pcBIdx]
        plot(x,y, xlab=pcCIdx, ylab=pcBIdx, cex.lab=0.8, cex.axis=0.8, bg=sampleColours, col=sampleLcolours, pch=samplePch, lwd=sampleLwd, cex=sampleSize)
        x <- sampleMeta[,pcAIdx]
        y <- sampleMeta[,pcCIdx]
        plot(x,y, xlab=pcAIdx, ylab=pcCIdx, cex.lab=0.8, cex.axis=0.8, bg=sampleColours, col=sampleLcolours, pch=samplePch, lwd=sampleLwd, cex=sampleSize)
        plot(x,y,xaxt="n",xlab="", ylab="",yaxt="n",type="n",bty="n")

        graphics::legend("topleft", ncol=1, inset=0.05, cex=1.0, 
                         legendData$label, pt.bg=legendData$colour, col=legendData$lcolour, 
                         pch=as.numeric(legendData$pch), pt.lwd=as.numeric(legendData$lwd), 
                         pt.cex=as.numeric(legendData$size))
    } else {
        #
        # Do a separate legend- useful if the names are long
        #
        if (!is.null(legendData)) {
            grDevices::dev.off()
            plotFilename  <- paste("pca-",plotName,"-legend", sep="")
            graphicFilenameRoot  <- paste(plotsFolder, plotFilename, sep="/")
            initializeGraphics (getGraphicsFilename (graphicFilenameRoot), params)
            graphics::plot.new()
            graphics::par(mar=c(0,0,0,0))
            graphics::legend("topleft", ncol=1, inset=0.05, cex=1.0, 
                         legendData$label, pt.bg=legendData$colour, col=legendData$lcolour, 
                         pch=as.numeric(legendData$pch), pt.lwd=as.numeric(legendData$lwd), 
                         pt.cex=as.numeric(legendData$size))
        }
    }
    grDevices::dev.off()
}

