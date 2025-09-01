
UI.SUPPORTED_MAP_TYPES <- c("location", "sampleCount", "drug", "mutation", "alleleProp", "diversity")

ui.makeSingleMap <- function (ctx, sampleSet, type,
                              #timePeriods=NULL,	TBD
                              measure=NULL,
                              aggregate="Province", minAggregateCount=1,
                              markerSize=c(4,40), markerColours="red3", markerFontSize=6,
                              colourBy="Province",
                              showNames=TRUE, nameFontSize=5,
                              ...) {

    if (!(type %in% UI.SUPPORTED_MAP_TYPES)) {
        stop (paste("Unsupported map type for UI:", type))
    }
    task <- paste0("map/", type)
    config <- context.getConfig(ctx)

    params <- param.makeParameterList (ctx, task, 
                  timePeriods=NULL, aggregate=aggregate, minAggregateCount=minAggregateCount, 
                  markerSize=markerSize, markerColours=markerColours, markerFontSize=markerFontSize, 
                  colourBy=colourBy, showNames=showNames, nameFontSize=nameFontSize, 
                  ...)											#; print(as.list(params))
    #
    # measure param
    #
    ui.checkSingleValue (measure, "measure")
    mArg <- list(measure=measure)
    if (type=="location") {
        params$analysis.measures <- "Location"
        
    } else if (type=="sampleCount") {
        params$analysis.measures <- "SampleCount"
        
    } else if (type=="drug") {
        params$analysis.measures <- param.getArgParameter (mArg, "measure", type="character", 
                                                           validValues=c("ALL", setup.getFeatureNames(config$drugPredictionFeatures)))

    } else if (type=="mutation") {
        params$analysis.measures <- param.getArgParameter (mArg, "measure", type="character", 
                                                           validValues=c("ALL", setup.getFeatureNames(config$drugMutationFeatures)))
        
    } else if (type=="alleleProp") {
        params$analysis.measures <- param.getArgParameter (mArg, "measure", type="character", 
                                                           validValues=c("ALL", setup.getFeatureNames(config$countableFeatures), 
                                                                         setup.getFeatureNames(config$amplificationFeatures), 
                                                                         setup.getFeatureNames(config$drugLocusFeatures)))
        
    } else if (type=="diversity") {
        params$analysis.measures <- param.getArgParameter (mArg, "measure", type="character", 
                                                           validValues=c("ALL", diversity.getDiversityMeasures()))
        
    }

    plotInfo <- execute.makeMapOnSampleSet (userCtx=ctx, sampleSetName=sampleSet, task=task, params=params)
    plotInfo
}
#
#
#
ui.checkSingleValue <- function (value, paramName) {
    if (length(value) > 1) {
        stop (paste(paramName, "can only accet a single value"))
    }
}
#
#
execute.makeMapOnSampleSet <- function (userCtx, sampleSetName, task, params) {
    #
    # Execute computations and plots 
    #
    execute.checkSampleSet (userCtx, sampleSetName)
    #
    # Execute the task
    #
    t <- execute.parseTask (task)
    if (t$name != "map") {
        stop (paste ("Invalid analysis task:", t$taskName))
    }
    plotInfo <- map.makeSingleMap (userCtx, sampleSetName, t$method, params)
    plotInfo
}
#
#
#
map.makeSingleMap <- function(userCtx, sampleSetName, mapType, params) {		#;print(mapType)	;print(params)
    #
    # Get the set of specs for the maps to be produced
    #
    mapSpecs <- map.createMapSpecs (userCtx, sampleSetName, mapType, params)
    mapSpec <- mapSpecs[[1]]
    #
    # Generate the ggplot2 plot object
    #
    mapPlot <- map.generateMapPlot (mapSpec, params)
    #
    # Add the "look and feel" elements to the map plot
    #
    mapPlot <- mapPlot + mapSpec$master$theme
    #
    # Deal with the legend - TBD
    #
    plotInfo <- map.processLegend (mapPlot, mapSpec, params)
    plotInfo
}
