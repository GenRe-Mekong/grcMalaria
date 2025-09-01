###################################################################
# Task Parameters retrieval
###################################################################
#
param.getParam <- function (paramName, params) {
    if (is.null(params)) {
        stop ("Missing parameters list")
    }
    value <- params[[paramName]]
    if (is.null(value)) {
        stop (paste("Missing parameter:",paramName))
    }
    value
}
#
###################################################################
# Argument parsing and translating to parameters
###################################################################
#
param.getArgParameter <- function (args, argName, type="character", defaultValue=NULL, multiValue=FALSE, validValues=NULL, scaleFactor=NULL) {
    value <- args[[argName]]				#; print("******************") ; print(argName)
    if (is.null(value)) {
        value <- defaultValue
    }							
    if (!is.null(value)) {
        storage.mode(value) <- type			#; print(type)
    }							#; print(value)
    if (!multiValue) {
        if (length(value) > 1) {			#; print(paste("Multi-value:",length(value)))
            stop(paste0("Multiple values given for single-value argument '", argName, "'"))
        }		
    }
    if (!is.null(validValues)) {			#; print(paste("Valid values:",paste(validValues, collapse=",")))
        isValid <- TRUE
        if (is.null(value)) {  
            isValid <- FALSE
        } else if (multiValue) {
            isValid <- all(value %in% validValues)
        } else {
            isValid <- value %in% validValues
        }
        if (!isValid) {
            stop(paste0("Invalid value for argument '", argName, "': ", args[[argName]]))
        }
    }
    if ((type == "numeric") && (!is.null(scaleFactor))) {
        value <- value / scaleFactor			#; print(paste("Scaled value:",value))
    }
    value
}
#
###################################################################
# Parameters list creation, with defaults
###################################################################
#
param.makeParameterList <- function (ctx, task, ...) {
    args <- list(...)					#; print(args)
    t <- execute.parseTask (task)			#; print(t)
    taskName <- t$name
    method   <- t$method
    #
    # Create the list of parameters
    #
    p <- new.env ()
    #
    # Poulate according to the task
    #
    if (taskName == "map") {
        param.addPlotGeomParameters (p, args, ...)
        param.addMapParameters (ctx, p, args, method, ...)
        
    } else if (taskName %in% c("pca", "tree")) {	#; print(args$plots)
        param.addPlotGeomParameters (p, args, ...)
        p$plot.plotList <- param.getArgParameter (args, "plots", type="list", defaultValue=NULL, multiValue=TRUE)
        p$phylo.impute  <- param.getArgParameter (args, "impute", type="logical", defaultValue=FALSE)
        
    } else if (taskName == "graph") {
        param.addPlotGeomParameters (p, args, ...)
        p$cluster.clusterSet.name <- param.getArgParameter (args, "clusterSet",  type="character")
        p$graph.layoutAlgorithm   <- param.getArgParameter (args, "graphLayout", type="character", defaultValue="fr", validValues=c("fr"))
        p$graph.weightPower       <- param.getArgParameter (args, "weightPower", type="numeric",   defaultValue=2.0)
        
    } else {
        stop(paste("Invalid analysis task:", taskName))
    }
    p
}
#
#
#
param.addPlotGeomParameters <- function (p, args, ...) {
    #
    # First get the plot geometry parameters, which are needed to scale some of the plotting parameters
    #
    p$plot.width  <- param.getArgParameter (args, "width",  type="numeric", defaultValue=15)
    p$plot.height <- param.getArgParameter (args, "height", type="numeric", defaultValue=15)
    p$plot.units  <- param.getArgParameter (args, "units",  defaultValue="in", validValues=c("in", "cm", "mm", "px"))
    p$plot.dpi    <- param.getArgParameter (args, "dpi",    type="numeric", defaultValue=300)
    p$plot.file.format <- param.getArgParameter (args, "format", defaultValue="png", validValues=c("png", "pdf", "jpg"))
    #
    # Determine the scale factor needed to match the size of the device
    #
    p$plot.aspectRatio <- p$plot.width / p$plot.height
    #
    # Get the size requested in inches and get the scaling ratio that will make the smallest side equal to 16 inches
    #
    units <- p$plot.units; width <- p$plot.width; height <- p$plot.height; dpi <- p$plot.dpi
    if        (units == "in") {  inW <- width;      inH <- height
    } else if (units == "cm") {  inW <- width/2.54; inH <- height/2.54
    } else if (units == "mm") {  inW <- width/25.4; inH <- height/25.4
    } else if (units == "px") {  inW <- width/dpi;  inH <- height/dpi
    }
    inMax <- min(inW, inH)
    p$plot.scaleFactor <- 16 / inMax								#; print(p$plot.scaleFactor)
    #
    # Legend positioning
    #
    p$plot.legend.pos   <- param.getArgParameter (args, "legendPosition", defaultValue="inset", validValues=c("inset", "separate"))
    p$plot.legend.dir   <- param.getArgParameter (args, "legendDirection", defaultValue="vertical", validValues=c("vertical", "horizontal"))
}
#
#
#
param.addMapParameters <- function (ctx, p, args, taskMethod, ...) {
    config <- context.getConfig(ctx)
    scale <- p$plot.scaleFactor
    #
    # Aggregation levels- apply to all map types
    #
    aggLabels <- param.getArgParameter (args, "aggregate", defaultValue="Province", multiValue=TRUE, validValues=c("Province","District"))
    p$aggregation.levels <- map.getAggregationLevelsFromLabels (aggLabels)			#; print(p$aggregation.levels)
    if (taskMethod == "location") {
        defaultAggregateCountMin <- 1
    } else {
        defaultAggregateCountMin <- 10
    }
    p$map.aggregateCountMin <- param.getArgParameter (args, "minAggregateCount", type="numeric", defaultValue=defaultAggregateCountMin)
    #
    # Get the time slice intervals
    #
    if (taskMethod %in% c("sampleCount", "drug", "mutation", "alleleProp", "diversity")) {
          timePeriods <- param.getArgParameter (args, "timePeriods", type="list", defaultValue=NULL, multiValue=TRUE)
          p$analysis.timeIntervals <- parseTimeIntervals(timePeriods)				#; print(p$analysis.timeIntervals)
    }
    #
    # Markers and geographical names
    #
    if (taskMethod %in% c("location", "sampleCount", "drug", "mutation", "alleleProp", "diversity", "connect", "clusterPrevalence")) {
        defMarkerSize <- ifelse (taskMethod == "connect", 6, 16)				#; print (paste("defMarkerSize:",defMarkerSize))
        if (taskMethod == "alleleProp") {
            # For pie chartts, do not scale
            p$map.markerSize <- param.getArgParameter (args, "markerSize", type="numeric", multiValue=TRUE, defaultValue=defMarkerSize)
        } else {
            p$map.markerSize <- param.getArgParameter (args, "markerSize", type="numeric", multiValue=TRUE, defaultValue=defMarkerSize, scaleFactor=scale)
        }											#; print(p$map.markerSize)
        if (!(taskMethod %in% c("location", "connect"))) {
            p$map.markerValueFontSize <- param.getArgParameter (args, "markerFontSize", type="numeric", defaultValue=6, scaleFactor=scale)
        }
    }
    p$map.markerNames <- param.getArgParameter (args, "showNames", type="logical", defaultValue=FALSE)
    p$map.admNameLabelFontSize <- param.getArgParameter (args, "nameFontSize", type="numeric", defaultValue=5, scaleFactor=scale)
    #
    # Other scalable plotting parameters (we're not advertising how to control these initially) 
    #
    p$themeBaseSize        <- param.getArgParameter (args, "themeBaseSize",        type="numeric", defaultValue=11,  scaleFactor=scale)
    p$plotTitleSize        <- param.getArgParameter (args, "plotTitleSize",        type="numeric", defaultValue=1.2, scaleFactor=scale)
    p$axisTitleSize        <- param.getArgParameter (args, "axisTitleSize",        type="numeric", defaultValue=1.0, scaleFactor=scale)
    p$baseBorderWidth      <- param.getArgParameter (args, "baseBorderWidth",      type="numeric", defaultValue=1.0, scaleFactor=scale)
    p$markerBorderWidth    <- param.getArgParameter (args, "markerBorderWidth",    type="numeric", defaultValue=1.5, scaleFactor=scale)
    p$admNameLabelPadding  <- param.getArgParameter (args, "admNameLabelPadding",  type="numeric", defaultValue=0.3, scaleFactor=scale)
    p$admNameLineWidth     <- param.getArgParameter (args, "admNameLineWidth",     type="numeric", defaultValue=1.0, scaleFactor=scale)
    p$pieLineWidth         <- param.getArgParameter (args, "pieLineWidth",         type="numeric", defaultValue=0.5, scaleFactor=scale)
    p$legendBorderWidth    <- param.getArgParameter (args, "legendBorderWidth",    type="numeric", defaultValue=0.5, scaleFactor=scale)
    p$legendFontSize       <- param.getArgParameter (args, "legendFontSize",       type="numeric", defaultValue=4,   scaleFactor=scale)
    p$legendKeySize        <- param.getArgParameter (args, "legendKeySize",        type="numeric", defaultValue=1.2, scaleFactor=scale)
    p$connectCurveWidthMin <- param.getArgParameter (args, "connectCurveWidthMin", type="numeric", defaultValue=0.25, scaleFactor=scale)
    p$connectCurveWidthMax <- param.getArgParameter (args, "connectCurveWidthMax", type="numeric", defaultValue=4.0, scaleFactor=scale)
    #
    # Non-scalable and non-numeric parameters plotting parameters (we're not advertising how to control these initially) 
    #
    p$plotTitleSize   <- param.getArgParameter (args, "plotTitleSize",   type="numeric", defaultValue=1.2)
    p$axisTitleSize   <- param.getArgParameter (args, "axisTitleSize",   type="numeric", defaultValue=1.0)
    p$adm0BorderWidth <- param.getArgParameter (args, "adm0BorderWidth", type="numeric", defaultValue=1.5)
    p$adm1BorderWidth <- param.getArgParameter (args, "adm1BorderWidth", type="numeric", defaultValue=1.0)
    p$riverWidth      <- param.getArgParameter (args, "riverWidth",      type="numeric", defaultValue=1.0)
    p$lakeshoreWidth  <- param.getArgParameter (args, "lakeshoreWidth",  type="numeric", defaultValue=1.0)
    #
    # Plot-specific parameters
    #
    if (taskMethod == "location") {
        p$analysis.measures <- "Location"
        colourByLabel <- param.getArgParameter (args, "colourBy",  defaultValue="Province", validValues=c("Country", "Province"))
        p$map.markerColourAggLevel <- map.getAggregationLevelsFromLabels (colourByLabel)

    } else if (taskMethod == "sampleCount") {
        p$analysis.measures <- "SampleCount"
        p$map.aggregateCountMin <- param.getArgParameter (args, "minAggregateCount", type="numeric", defaultValue=1)
        colourByLabel <- param.getArgParameter (args, "colourBy",  defaultValue="Province", validValues=c("Country", "Province"))
        p$map.markerColourAggLevel <- map.getAggregationLevelsFromLabels (colourByLabel)

    } else if (taskMethod == "drug") {
        p$analysis.measures <- param.getArgParameter (args, "drugs", type="character", multiValue=TRUE, defaultValue="ALL", 
                                                      validValues=c("ALL",setup.getFeatureNames(config$drugPredictionFeatures)))
                                                      
    } else if (taskMethod == "mutation") {
        p$analysis.measures <- param.getArgParameter (args, "mutations", type="character", multiValue=TRUE, defaultValue="ALL", 
                                                      validValues=c("ALL",setup.getFeatureNames(config$drugMutationFeatures)))
        p$map.markerColours <- param.getArgParameter (args, "markerColours", multiValue=TRUE, defaultValue="red3")

    } else if (taskMethod == "alleleProp") {
        p$analysis.measures <- param.getArgParameter (args, "mutations", type="character", multiValue=TRUE, defaultValue="ALL", 
                                                      validValues=c("ALL", setup.getFeatureNames(config$countableFeatures), 
                                                                     setup.getFeatureNames(config$amplificationFeatures), 
                                                                     setup.getFeatureNames(config$drugLocusFeatures)))
        p$map.alleleColours <- param.getArgParameter (args, "alleleColours", multiValue=TRUE, defaultValue=NULL)

    } else if (taskMethod == "diversity") {
        p$analysis.measures <- param.getArgParameter (args, "measures", type="character", multiValue=TRUE, defaultValue="ALL", 
                                                      validValues=c("ALL",diversity.getDiversityMeasures()))
        p$map.markerColours <- param.getArgParameter (args, "markerColours", multiValue=TRUE, defaultValue="red3")
        
    } else if (taskMethod == "connect") {
        p$analysis.measures <- param.getArgParameter (args, "measures", type="character", multiValue=TRUE, defaultValue="ALL", 
                                                      validValues=c("ALL",connectMap.getConnectednessMeasures()))
        p$map.connect.identity.min     <- param.getArgParameter (args, "minIdentity", type="numeric", multiValue=TRUE, defaultValue=1.0)
        p$map.connect.meanDistance.min <- param.getArgParameter (args, "meanDistanceLevels", type="numeric", multiValue=TRUE, defaultValue=0.5)

    } else if (taskMethod == "barcodeFrequency") {
        p$map.cluster.visualizations <- param.getArgParameter (args, "type", type="character", multiValue=TRUE, defaultValue=c("bar","pie"), 
                                                               validValues=c("bar","pie"))
        p$map.cluster.markerScale <- param.getArgParameter (args, "markerScale", type="numeric", defaultValue=0.8)
        p$map.cluster.markerSampleCount <- param.getArgParameter (args, "markerSampleCount", defaultValue="mean")

    } else if (taskMethod == "clusterSharing") {
        p$cluster.clusterSet.name <- param.getArgParameter (args, "clusterSet", type="character")
        p$map.cluster.visualizations <- param.getArgParameter (args, "type", type="character", multiValue=TRUE, defaultValue=c("bar","pie"), 
                                                               validValues=c("bar","pie"))
        p$map.cluster.markerScale <- param.getArgParameter (args, "markerScale", type="numeric", defaultValue=0.8)
        p$map.cluster.markerSampleCount <- param.getArgParameter (args, "markerSampleCount", defaultValue="mean")

    } else if (taskMethod == "clusterPrevalence") {
        p$cluster.clusterSet.name <- param.getArgParameter (args, "clusterSet", type="character")

    } else {
        stop(paste("Invalid plot task:", paste0("map/",taskMethod)))
    }
}

