#############################################################
#
#' grcMalaria: Genetic Epidemiology Analysis using SpotMalaria Genetic Report Cards
#'
#'\if{html}{\figure{logo.jpg}{options: align='right' alt='logo' width='120'}}
#' A user-friendly, open-source R package, designed to make genetic epidemiology analysis tasks accessible.
#' This package facilitates the translation of genetic information derived from malaria parasites from SpotMalaria Genetic Report Cards (GRC),
#' into intuitive geographical maps of prevalence, diversity, relatedness.
#' This software library is also capable of identifying circulating strains, characterizing drug resistance profiles, and mapping spread.
#' The R documentation was updated to align with V2.0.0
#' Refer to the user-guide (link below) for the latest instructions on how to use this package.
#'
#' @section Author:
#' Olivo Miotto \email{olivo@tropmedres.ac}
#'
#' @seealso Useful links:
#' \itemize{
#'  \item user guide {\url{https://genremekong.org/tools/grcmalaria-r-package-user-guide}}
#' }
#'
#' @name grcMalaria
"_PACKAGE"

#############################################################
#
#' Load a Genetic Report Card data file
#'
#' Load a GRC data file (Microsoft Excel Worksheet format), by indicating the file path in the loadGrc() function.
#' Save the GRC onto your computer in order to load onto grcMalaria
#' All types of GRC data files are supported, whether it is a public data file, or privately owned GRC data for study owners.
#' A GRC data file consists of groups of variables, including metadata fields with geographical and spatial data (grey headers),
#' private metadata (pink, not included in public data), predicted drug resistance (orange),
#' sample characteristics including species determination, and complexity of infection (light blue),
#' drug resistance amino acid haplotypes (green), amino acid variant for individual drug resistance variants (blue),
#' quantitative polymerase chain reaction (qPCR) results for gene copy number variation (CNV) for plasmepsin2/3 (pm23),
#' and multidrug resistance protein 1 (Pfmdr1) (purple), and lastly genotypes used for genetic barcodes (red).
#'
#' Note that, from V2.0 onwards, GRC files contain a Properties sheet, which describes the contents of the other sheets.
#' Therefore, parameters species, version and sheet are ignored for V2.0 or greater. Those paramenters are required for GRC V1.x.
#' 
#' @param file The path to the GRC data Microsoft Excel file (use forward slashes in the path)
#' @param species The species being analysed (only valid for GRC versions 1.x, can be "Pf" or "Pv").
#' @param version The version number of the GRC data file format (only valid for GRC versions 1.x)
#' @param sheet Name of worksheet in the Excel file that contains the GRC data (only valid for GRC versions 1.x; from V2.0 onwards, the data sheet can only be named "Data")
#'
#' @return A list containing a data frame with the data ready to be analyzed, plus some configuration metadata
#' @export
#'
#' @examples \dontrun{
#' # Load data file
#' # Change the path to where your file is located before running the code
#' Data <- loadGrc("D:/myFiles/name_file.xlsx",
#'                sheet = "GRC", species = "Pf", version = "1.4")
#'}
#'
#
loadGrc <- function (file, species=NULL, version=NULL, sheet=NULL) {
    grc <- grcData.load (file, species, version, sheet)
    grc
}

#############################################################
#
#' Combine two Genetic Report Card (GRC) datasets
#'
#' Merge two GRC datasets which have been loaded from data files. The new dataset is added to the source dataset.
#' If the new dataset contains new columns, they can be either ignored, or the same columns are added to the source dataset,
#' with null values. If the new dataset contains sample data that is already in the source dataset, they can be either ignored,
#' or the source data rows can be overwritten.
#'
#' @param srcGrc The source GRC dataset (loaded by loadGrc())
#' @param newGrc The new GRC dataset to be added(loaded by loadGrc())
#' @param overwrite Boolean: if TRUE, sample rows in newGrc will replace rows in srcGrc that have the same sampleId. If FALSE, they will be disregarded.
#' @param extendColumns Boolean: if TRUE, new columns in newGrc will also be added to srcGrc, filled with <NA> values.  If FALSE, they will be disregarded.
#'
#' @return A list containing a data frame with the merged data ready to be analyzed, plus some configuration metadata
#' @export
#'
#' @examples \dontrun{
#' ## Load dataset1 ##
#' Sheet1 <- loadGrc("D:/.../sheet1.xlsx",
#'                  sheet="GRC", species="Pf", version="1.4")
#' ## Load dataset2 ##
#' Sheet2 <- loadGrc("D:/.../sheet2.xlsx",
#'                  sheet="GRC", species="Pf", version="1.4")
#' ## Merge datasets ##
#' mergeData <- mergeGrc(Sheet1, Sheet2)
#'}
#
mergeGrc <- function (srcGrc, newGrc, overwrite=FALSE, extendColumns=FALSE) {
  grc <- grcData.merge (srcGrc, newGrc, overwrite, extendColumns)
  grc
}

#############################################################
#
#' Get a list of items that con be measured and mapped
#'
#' The sample features you can analyze depends on the format of each specific GRC. This function can be used
#' to discover what features can be analyzed for the GRC you have loaded. For example, you can ask for a 
#' list of "mutation" features, which will return a list of names of mutations that are
#' supported by the GRC. Then you can supply one or more of these names to the "mutations" parameter of
#' function mapMutationPrevalence(), to botain prevalence maps.
#'
#' @param grc The GRC to provides the data to be analyzed (loaded by loadGrc())
#' @param category String: The category of features to be listed. Can be one of ["drug", "locus", "mutation", "amplification"]
#'
#' @return A vector of strings, consisting of the names of the  features of the specified category that are supported by this GRC.
#' @export
#'
#' @examples \dontrun{
#' ## Load dataset ##
#' grc <- loadGrc("D:/.../sheet1.xlsx", species="Pf", version="2.0")
#' mutNames <- getGrcFeatureNames (grc, "mutation")
#' ... 
#' mapMutationPrevalence(ctx, sampleSet="set1", mutations=mutNames[1:3], aggregate="District")
#'}
#
getGrcFeatureNames <- function (grc, category) {
  featureNames <- grcData.getFeatureNames (grc, category)
  featureNames
}

#############################################################
#
#' Initialize a Genetics Report Card (GRC) Dataset
#'
#' This function creates a datafolder for subsequent analysis, removes samples and SNPs with high missingness, based on the thresholds specified.
#' It then performs barcode imputation and calculates genetic distance matrices using the remaining data.
#' Filtering is performed to minimise errors that can occur due to insufficient DNA concentrations for sequencing.
#' Insufficiency in DNA could be attributed to low parasitemia, or incorrect Plasmodium species classification.
#' By default, SNPs with missing calls in more than 20% of samples (minSnpTypability=0.8),
#' and samples with missing alleles in more than 25% of the remaining SNPs (minSampleTypability=0.75), are removed.
#'
#' The function creates two directories, ‘data’, and ‘out’, within the specified path.
#' The ‘data’ directory contains a subdirectory for each GRC dataset. Within this subdirectories, four folders are created:
#' (i) ‘barcode’: contains filtered and imputed barcodes in tabular and fasta format
#' (ii) ‘distance’: contains a filtered and an imputed distance matrix
#' (iii) ‘genotype’: contains filtered and imputed genotypes
#' (iv) ‘meta’: contains unfiltered, filtered, and imputed sample metadata
#'
#' @param grc The data obtained from reading the GRC data with loadGrc().
#' @param dir You will want to change the path and folder name to where you want the output file to be stored on your computer.
#' @param outDir The subfolder of dir where the analysis outputs (tables and maps) will be written. Default: subfolder "out"
#' @param minSnpTypability The minimum proportion of non-missing samples for a barcode position to be retained in analysis. The default is 0.8. This means SNPs where the allele is missing in more than 20 percent of samples are removed.
#' @param minSampleTypability he minimum proportion of non-missing positions for a sample to be retained in analysis. The default is 0.75. This means samples where more than 25 percent of the remaining SNPs are missing are removed.
#' @param maxImputedProportion The maximum proportion of barcode positions to be imputed for a filtered sample to be included in the imputed dataset.
#' @param clearCacheData Removes previously stored data from analyzing the same dataset.
#'
#' @return A analysis context object, which is a list that contains all the data for analysis, which will be passed to subsequent analysis tasks.
#' @export
#'
#' @examples \dontrun{
#' # Change the path to where you want output file to be
#' ctx <- initializeContext(Data,
#'                               dir="D:/...", #Change the path to where you want output file to be
#'                               minSnpTypability=0.8, minSampleTypability=0.75)
#'}
#
initializeContext <- function (grc, dir=".", outDir="out",
                               minSnpTypability=0.8,
                               minSampleTypability=0.75,
                               maxImputedProportion=(1-minSampleTypability),
                               clearCacheData=FALSE) {
  options(scipen=10)
  options(stringsAsFactors=FALSE)

  config <- setup.getConfig (grc, dir, outDir, minSnpTypability, minSampleTypability, maxImputedProportion)
  ctx <- context.createRootContext (grc$data, config, clearCacheData)
  ctx
}

#############################################################
#
#' Select a set of samples using metadata
#'
#' Selects a set of sample for analysis, based on their metadata values.
#' A given analysis context can contain multiple sample sets. These can be labelled with different names.
#' The selection criteria are specified as a list containing a sequence of lists with two elements each:
#' "field" which is the column name to be checked, and "values" which is an array of possible values that can be matched for a sample to be selected.
#' A selected sample must match all the criteria.
#' Repeat this step for the different sample sets you would like to analyze
#'
#' @param ctx The analysis context, created by intializeContext().
#' @param sampleSetName Select a name for the sample set. You will use this for subsequent analysis tasks. This will create a folder with the same name e.g. “Laos”.
#' @param select The criteria for sample selection. Selected samples will match all specified criteria. The values must be a list. "field" corresponds to a column name in the datafile. "values" correspond to an array of possible values that can be matched for a sample to be selected. Use "," to add additional parameters inside a quotation mark " ".
#'
#' @return The analysis context, augmented with the new sample set
#' @export
#'
#' @examples \dontrun{
#' ## Select a sample set to work on ##
#'
#' # To select samples from 1 field (1 column in the GRC)
#' selectSampleSet(ctx, sampleSetName="EBKK", select=list(
#'                                 list(field="Country", values=c("VN", "KH", "LA")) ))
#'
#' # To select samples from 2 fields
#' selectSampleSet(ctx, sampleSetName="Laos", select=list(
#'                                 list(field="TimePoint", values=c("D00H00","-")),
#'                                 list(field="Study", values=c("1208-PF-LA-CMPE-GENRE")) ))
#'
#' # To select samples from 3 fields
#' selectSampleSet(ctx, sampleSetName="SouthLA_2017", select=list(
#'                                 list(field="Country", values="LA"),
#'                                 list(field="AdmDiv1", values=c("Attapeu", "Champasak")),
#'                                 list(field="Year", values=c("2017", "2018")) ))
#'
#'}
#
selectSampleSet <- function (ctx, sampleSetName, select) {
  sampleSet.selectSampleSet (ctx, sampleSetName, select)
}
#############################################################
#
#' Map of Locations
#'
#' Creates maps showing markers for all Provinces/Districts/Sites where samples were collected.
#' For each location, we place a marker on the map, coloured according to the administrative division it is in.
#'
#' @param ctx The analysis context, created by intializeContext().
#' @param sampleSet The name of the sample set being used, which must have been previously created by selectSampleSet().
#' @param aggregate The administrative level at which we aggregate (Province/District/Site). Separate maps are created for each administrative level.
#' @param markerSize Allows adjustment of the size of markers on the map. If only one value is passed, all markers will be of that size; if two values are passed, they will be used as the min and max size of the marker, whose size will reflect the number of samples.
#' @param colourBy Shows the aggregation level to be used to colour the markers (”Province” or “Country”). For each aggregation unit, we place a marker on the map, coloured according to the level of resistance to the drug or mutation, with a label indicating the prevalence.
#' @param showNames If TRUE, labels are shown with the name of the aggregation unit (Province or District)
#' @param nameFontSize Set the font size of the geographical name labels outside the markers on the map, if showNames=TRUE, default=5
#' @param ... Aesthetics: any of the following plot parameters,
#'            width: the width of the plot (numeric), default=15
#'            height: the height of the plot (numeric), default=15
#'            units: the units in which the width and height are expressed. Supported values: "in" (default), "cm", "mm", "px"
#'            dpi: the resolution of the plot output, expressed as dots per inch, default=300
#'            format: the file format in which the plot will be saved. Supported values: "png" (default), "pdf"
#'            legendPosition: specifies where the legend should be plotted. Supported values: "inset" (default), "separate"
#'            legendWidth: specifies how wide a fixed width legend space should be, default="NULL"
#'            legendDirection: specifies location of legend. Supported values: "vertical" (default), "horizontal"
#'            legendFontSize: specifies font size of the legend (numeric), default=4
#'            axisTitleSize: specifies axis label font size (numeric), default=1#'
#' @export
#'
#' @examples \dontrun{
#' # Create 2 maps showing the locations for sampleSet Cambodia on a Province, and a District level
#' mapLocations (ctx, sampleSet="Cambodia", aggregate=c("Province", "District"),
#'              markerSize=5, colourBy="Country", showNames=FALSE,
#'              width=18, height=10, units="in", format="pdf")
#'}
#
mapLocations <- function (ctx, sampleSet,
                          aggregate="Province",
                          markerSize=c(4,40), colourBy="Province",
                          showNames=TRUE, nameFontSize=5,
                          ...) {

  task <- "map/location"
  params <- param.makeParameterList (ctx, task,
                                     aggregate=aggregate, markerSize=markerSize, colourBy=colourBy, showNames=showNames, nameFontSize=nameFontSize,
                                     ...)
  execute.executeOnSampleSet (userCtx=ctx, sampleSetName=sampleSet, task=task, params=params)
}
#############################################################
#
#' Map of Sample Counts
#'
#' Creates maps showing the numbers of samples collected at different aggregation levels.
#' Number of samples are shown on the marker.
#' If number of samples is less than minimum count of aggregated samples (minAggregateCount), the marker will not appear.
#' File name with 'filtered’ is showing samples that passed the quality filtering threshold (default is minSampleTypability=0.75, meaning samples that have more than 25% missingness in their barcode are filtered out and not shown on the map).
#' File name with ‘unfiltered’ is showing all the samples without applying quality filtering.
#'
#' @param ctx The analysis context, created by intializeContext().
#' @param sampleSet The name of the sample set being used, which must have been previously created by selectSampleSet().
#' @param timePeriods Time-sequence maps can be implemented using a parameter called timePeriods parameter. When this is passed to a function, it will partition samples into time-interval plots. Time intervals parameters are available in 5 analyses: mapSampleCounts(), mapDrugResistancePrevalence(), mapMutationPrevalence(), mapAlleleProportions() and mapDiversity().
#'                    Why use time-intervals parameter:
#'                    It will produce maps with consistent geographical boundaries
#'                    Users can slice up the dataset in any specified time period
#' @param aggregate The administrative level at which we aggregate (Province or District). Separate maps are created for each administrative level.
#' @param minAggregateCount The minimum count of aggregated samples. To avoid estimating on very small samples, one can set a minimum count of samples, below which the marker is not shown.
#' @param markerSize Allows adjustment of the size of markers on the map. If only one value is passed, all markers will be of that size; if two values are passed, they will be used as the min and max size of the marker, whose size will reflect the number of samples.
#' @param markerFontSize Allows adjustment of the font size shown on the markers (numeric), default=6
#' @param colourBy Shows the aggregation level to be used to colour the markers (”Province” or “Country”). For each aggregation unit, we place a marker on the map, coloured according to the level of resistance to the drug or mutation, with a label indicating the prevalence.
#' @param showNames If TRUE, labels are shown with the name of the aggregation unit (Province or District)
#' @param nameFontSize Allows adjustment of the Province or District label name font size (numeric), default=5
#' @param ... Aesthetics: any of the following plot parameters,
#'            width: the width of the plot (numeric), default=15
#'            height: the height of the plot (numeric), default=15
#'            units: the units in which the width and height are expressed. Supported values: "in" (default), "cm", "mm", "px"
#'            dpi: the resolution of the plot output, expressed as dots per inch, default=300
#'            format: the file format in which the plot will be saved. Supported values: "png" (default), "pdf"
#'            legendPosition: specifies where the legend should be plotted. Supported values: "inset" (default), "separate"
#'            legendWidth: specifies how wide a fixed width legend space should be, default="NULL"
#'            legendDirection: specifies location of legend. Supported values: "vertical" (default), "horizontal"
#'            legendFontSize: specifies font size of the legend (numeric), default=4
#'            axisTitleSize: specifies axis label font size (numeric), default=1
#'
#' @export
#'
#' @examples \dontrun{
#' ## To implement time-sequence maps ##
#' # 1. Specify time periods
#' # name: the name parameter is required and is added to the end of the file name of the produced 
#' # files.
#' # There are two types of time interval: type="year" uses start date and by default show until 
#' # the end of the defined year; type="period" allows user to specify any time period using both 
#' # the start and end date.
#' # start and end - start and end date. The parameter start must be a date in the 
#' # format dd-MMM-yyyy.
#' # NOTE: Public GRC data files contains only year of collection, therefore only the year in the 
#' # provided “start” date is used to produce time-sequence maps.
#' periods <- list(
#'                list(name="2021", type="year", start="1-Jan-2021"),
#'                list(name="2020", type="year", start="1-Jan-2020"),
#'                list(name="2017-19", type="period", start="1-Jan-2017", end="31-Dec-2019")
#'             )
#'
#'
#' # 2. Apply timePeriods parameter to mapSampleCounts
#' mapSampleCounts(ctx, sampleSet="Laos", timePeriods = periods,
#'                 aggregate=c("Province","District"), minAggregateCount=1,
#'                 markerSize=c(10,40), colourBy="Province", showNames=TRUE,
#'                 width=15, height=15)
#'}
#
mapSampleCounts <- function (ctx, sampleSet, timePeriods=NULL,
                             aggregate="Province", minAggregateCount=1,
                             markerSize=c(4,40), markerFontSize=6, colourBy="Province",
                             showNames=TRUE, nameFontSize=5,
                             ...) {

  task <- "map/sampleCount"
  params <- param.makeParameterList (ctx, task,
                                     timePeriods=timePeriods, aggregate=aggregate, minAggregateCount=minAggregateCount, markerSize=markerSize,
                                     markerFontSize=markerFontSize, colourBy=colourBy, showNames=showNames, nameFontSize=nameFontSize,
                                     ...)
  execute.executeOnSampleSet (userCtx=ctx, sampleSetName=sampleSet, task=task, params=params)
}
#
#############################################################
#
#' Map prevalence of Drug resistance (Predicted Phenotype)
#'
#' Creates a map showing the levels of resistance to a particular antimalarial drug for different administrative divisions.
#' Based on published genetic markers, each sample has a predicted phenotype for different types of antimalarial drugs.
#' If multiple drugs are specified, then different maps will be created for different drugs.
#'
#' @param ctx The analysis context, created by intializeContext()
#' @param sampleSet The name of the sample set being used, which must have been previously created by selectSampleSet()
#' @param timePeriods Time-sequence maps can be implemented using a parameter called timePeriods parameter. When this is passed to a function, it will partition samples into time-interval plots. Time intervals parameters are available in 5 analyses: mapSampleCounts(), mapDrugResistancePrevalence(), mapMutationPrevalence(), mapAlleleProportions() and mapDiversity().
#'                    Why use time-intervals parameter:
#'                    It will produce maps with consistent geographical boundaries
#'                    Users can slice up the dataset in any specified time period
#' @param drugs The antimalarial drugs for which prevalence of phenotypic resistance will be estimated; "ALL" creates maps for all the drugs for which phenotypic resistance predictions are available,
#'              which include "Artemisinin", "Chloroquine", "Piperaquine", "DHA-PPQ" (i.e. Dihydroartemisinin/piperaquine), "Mefloquine", "Sulfadoxine", "Pyrimethamine", "S-P" (i.e. Sulphadoxine-Pyrimethamine),
#'              "AS-MQ" (i.e. Artesunate-Melfoquine), and "S-P-IPTp" (i.e. intermittent preventive treatment of malaria in pregnancy using Sulfadoxine-Pyrimethamine)
#'              To specify a drug put the drug name in between quotation marks e.g. "Artemisinin", or  c("Artemisinin", "Chloroquine", "S-P") to select several specific drugs.
#' @param aggregate The administrative level at which we aggregate. Separate maps are created for each administrative level.
#' @param minAggregateCount The minimum count of aggregated samples. To avoid estimating on very small samples, one can set a minimum count of samples, below which the marker is not shown.
#' @param markerSize Allows adjustment of the size of markers on the map. If only one value is passed, all markers will be of that size; if two values are passed, they will be used as the min and max size of the marker, whose size will reflect the number of samples.
#' @param markerFontSize Allows adjustment of the font size shown on the markers (numeric), default=6.
#' @param showNames If TRUE, labels are shown with the name of the aggregation unit (Province or District).
#' @param nameFontSize Allows adjustment of the Province or District label name font size (numeric), default=5.
#' @param ... Aesthetics: any of the following plot parameters,
#'            width: the width of the plot (numeric), default=15
#'            height: the height of the plot (numeric), default=15
#'            units: the units in which the width and height are expressed. Supported values: "in" (default), "cm", "mm", "px"
#'            dpi: the resolution of the plot output, expressed as dots per inch, default=300
#'            format: the file format in which the plot will be saved. Supported values: "png" (default), "pdf"
#'            legendPosition: specifies where the legend should be plotted. Supported values: "inset" (default), "separate"
#'            legendWidth: specifies how wide a fixed width legend space should be, default="NULL"
#'            legendDirection: specifies location of legend. Supported values: "vertical" (default), "horizontal"
#'            legendFontSize: specifies font size of the legend (numeric), default=4
#'            axisTitleSize: specifies axis label font size (numeric), default=1
#'
#' @export
#'
#' @examples \dontrun{
#' ## To implement time-sequence maps ##
#' # 1. Specify time periods
#' # name: the name parameter is required and is added to the end of the file name of the produced 
#' # files.
#' # There are two types of time interval: type="year" uses start date and by default show until 
#' # the end of the defined year; type="period" allows user to specify any time period using both 
#' # the start and end date.
#' # start and end - start and end date. The parameter start must be a date in the 
#' # format dd-MMM-yyyy.
#' # NOTE: Public GRC data files contains only year of collection, therefore only the year in the 
#' # provided “start” date is used to produce time-sequence maps.
#' periods <- list(
#'                list(name="2021", type="year", start="1-Jan-2021"),
#'                list(name="2020", type="year", start="1-Jan-2020"),
#'                list(name="2017-19", type="period", start="1-Jan-2017", end="31-Dec-2019")
#'             )
#'
#' # 2. Apply timePeriods parameter to mapDrugResistancePrevalence
#' mapDrugResistancePrevalence(ctx, sampleSet="Laos", timePeriods = periods,
#'                             drugs="ALL", aggregate=c("Province","District"),
#'                             minAggregateCount=10, showNames=TRUE, markerSize=16,
#'                             width=15, height=15)
#'
#' # How to read the output:
#' # Number on the marker shows the proportion of samples with resistant markers for the specified 
#' # drug, ranging between 0 and 1 where 0 means no samples carry the resistant marker and 1 means 
#' # 100% of the samples carry the markers used to predict resistance.
#' # Resistance prevalence is calculated by dividing the number of resistant samples by the sum of 
#' # resistant and sensitive samples, fr = Resistant samples/(Resistant + Sensitive sample)
#'}
#'
#' @seealso Useful links:
#' \itemize{
#' \item Details of how drug resistance phenotypes are predicted based on genetic markers {\url{http://ngs.sanger.ac.uk/production/malaria/Resource/29/20200705-GenRe-05-PhenotypeRules-0.39.pdf}}
#' \item Other useful background information about the SpotMalaria platform  {\url{https://www.malariagen.net/resource/29}}
#' }
#
mapDrugResistancePrevalence <- function (ctx, sampleSet, timePeriods=NULL,
                                         drugs="ALL",
                                         aggregate="Province", minAggregateCount=10,
                                         markerSize=16, markerFontSize=6,
                                         showNames=TRUE, nameFontSize=5,
                                         ...) {

  task <- "map/drug"
  params <- param.makeParameterList (ctx, task,
                                     timePeriods=timePeriods, drugs=drugs, aggregate=aggregate, minAggregateCount=minAggregateCount, markerSize=markerSize,
                                     markerFontSize=markerFontSize, showNames=showNames, nameFontSize=nameFontSize,
                                     ...)
  execute.executeOnSampleSet (userCtx=ctx, sampleSetName=sampleSet, task=task, params=params)
}

#############################################################
#
#' Map the prevalence of genetic mutations related to drug resistance
#'
#' This creates maps of frequency of variants known to be associated with drug resistance.
#' Mutation frequency is the ratio of known mutant variants divided by total variants in the population.
#'
#' @param ctx The analysis context, created by intializeContext().
#' @param sampleSet The name of the sample set being used; must have been previously created by selectSampleSet().
#' @param timePeriods Time-sequence maps can be implemented using a parameter called timePeriods parameter. When this is passed to a function, it will partition samples into time-interval plots. Time intervals parameters are available in 5 analyses: mapSampleCounts(), mapDrugResistancePrevalence(), mapMutationPrevalence(), mapAlleleProportions() and mapDiversity().
#'                    Why use time-intervals parameter:
#'                    It will produce maps with consistent geographical boundaries
#'                    Users can slice up the dataset in any specified time period
#' @param mutations The genetic mutation types for which prevalence will be calculated; "ALL" creates maps for all 42 available genetic mutation types.
#'                  Available mutation types include: "crt_C72S",	"crt_M74I",	"crt_N75E",	"crt_N75D",	"crt_K76T",	"crt_T93S",	"crt_H97Y",	"crt_H97L",	"crt_I218F",	"crt_A220S",
#'                  "crt_Q271E",	"crt_N326S",	"crt_N326D",	"crt_T333S",	"crt_I356T",	"crt_I356L",	"crt_R371I",	"dhfr_N51I",	"dhfr_C59R",	"dhfr_C59H",	"dhfr_S108N",
#'                  "dhfr_S108T",	"dhfr_I164L",	"dhps_S436A",	"dhps_S436F",	"dhps_A437G",	"dhps_K540E",	"dhps_K540N",	"dhps_A581G",	"dhps_A613T",	"dhps_A613S",	"mdr1_N86Y",
#'                  "mdr1_Y184F",	"mdr1_S1034I",	"mdr1_F1226Y",	"mdr1_D1246Y",	"arps10_V127M",	"arps10_D128H",	"arps10_D128Y",	"fd_D193Y",	"mdr2_T484I",	"exo_E415G".
#'                  To specify a mutation put the mutation type in between quotation marks e.g. "crt_N75E", or  c("arps10_D128Y","fd_D193Y","crt_I218F") to select several genetic mutation types.
#' @param aggregate The administrative level at which we aggregate. Separate maps are created for each administrative level.
#' @param minAggregateCount The minimum count of aggregated samples. To avoid estimating on very small samples, one can set a minimum count of samples, below which the marker is not shown.
#' @param markerSize Allows adjustment of the size of markers on the map. If only one value is passed, all markers will be of that size; if two values are passed, they will be used as the min and max size of the marker, whose size will reflect the number of samples.
#' @param markerFontSize Allows adjustment of the font size shown on the markers (numeric), default=6.
#' @param showNames If TRUE, labels are shown with the name of the aggregation unit (Province or District).
#' @param nameFontSize Allows adjustment of the Province or District label name font size (numeric), default=5.
#' @param ... Aesthetics: any of the following plot parameters,
#'            width: the width of the plot (numeric), default=15
#'            height: the height of the plot (numeric), default=15
#'            units: the units in which the width and height are expressed. Supported values: "in" (default), "cm", "mm", "px"
#'            dpi: the resolution of the plot output, expressed as dots per inch, default=300
#'            format: the file format in which the plot will be saved. Supported values: "png" (default), "pdf"
#'            legendPosition: specifies where the legend should be plotted. Supported values: "inset" (default), "separate"
#'            legendWidth: specifies how wide a fixed width legend space should be, default="NULL"
#'            legendDirection: specifies location of legend. Supported values: "vertical" (default), "horizontal"
#'            legendFontSize: specifies font size of the legend (numeric), default=4
#'            axisTitleSize: specifies axis label font size (numeric), default=1
#'
#' @export
#'
#' @seealso Useful links:
#' \itemize{
#'  \item user guide for information on functionality of specific genetic mutation types {\url{https://genremekong.org/tools/grcmalaria-r-package-user-guide}}
#' }
#'
#' @examples \dontrun{
#' ## To implement time-sequence maps ##
#' # 1. Specify time periods
#' # name: the name parameter is required and is added to the end of the file name of the produced 
#' # files.
#' # There are two types of time interval: type="year" uses start date and by default show until 
#' # the end of the defined year; type="period" allows user to specify any time period using both 
#' # the start and end date.
#' # start and end - start and end date. The parameter start must be a date in the 
#' # format dd-MMM-yyyy.
#' # NOTE: Public GRC data files contains only year of collection, therefore only the year in the 
#' # provided “start” date is used to produce time-sequence maps.
#' periods <- list(
#'                list(name="2021", type="year", start="1-Jan-2021"),
#'                list(name="2020", type="year", start="1-Jan-2020"),
#'                list(name="2017-19", type="period", start="1-Jan-2017", end="31-Dec-2019")
#'             )
#'
#' # 2. Apply timePeriods parameter to mapMutationPrevalence
#' mapMutationPrevalence(ctx, sampleSet="Laos", timePeriods = periods,
#'                       mutations="ALL", aggregate="District",
#'                       minAggregateCount=10, showNames=TRUE, markerSize=16,
#'                       width=15, height=15)
#'
#' # How to read the output:
#' # Number on the marker shows frequency of specified variant. The scale ranges between 0 and 1
#' #where 0 means no samples carry the specified allele and 1 means 100% of the samples carry the
#' #specified allele.
#' # Mutation frequency = selected variant / total variants in the population.
#'}
#
mapMutationPrevalence <- function (ctx, sampleSet, timePeriods=NULL,
                                   mutations="ALL",
                                   aggregate="Province", minAggregateCount=10,
                                   markerSize=16, markerFontSize=6,
                                   showNames=TRUE, nameFontSize=5,
                                   ...) {

  task <- "map/mutation"
  params <- param.makeParameterList (ctx, task,
                                     timePeriods=timePeriods, mutations=mutations, aggregate=aggregate, minAggregateCount=minAggregateCount,
                                     markerSize=markerSize, markerFontSize=markerFontSize, showNames=showNames, nameFontSize=nameFontSize,
                                     ...)
  execute.executeOnSampleSet (userCtx=ctx, sampleSetName=sampleSet, task=task, params=params)
}

#############################################################
#
#' Map the proportion of alleles for a given mutation (including Pfkelch13 and gene amplifications).
#'
#' This function creates maps with pie-charts showing the proportion of selected alleles in each administrative division.
#' Markers are displayed as pie-chart showing proportion of all the homozygous alleles for the specified gene.
#' Missing genotypes and heterozygous genotypes are excluded in the calculations.
#'
#' @param ctx The analysis context, created by intializeContext()
#' @param sampleSet The name of the sample set being used; must have been previously created by selectSampleSet()
#' @param timePeriods Time-sequence maps can be implemented using a parameter called timePeriods parameter. When this is passed to a function, it will partition samples into time-interval plots. Time intervals parameters are available in 5 analyses: mapSampleCounts(), mapDrugResistancePrevalence(), mapMutationPrevalence(), mapAlleleProportions() and mapDiversity().
#'                    Why use time-intervals parameter:
#'                    It will produce maps with consistent geographical boundaries
#'                    Users can slice up the dataset in any specified time period
#' @param mutations The genetic mutation types for which prevalence will be calculated; "ALL" creates maps for all 11 available genetic mutations
#'                  Available genetic mutations include: "Pfkelch13", "pm23-Amp", "mdr1-Amp", "crt", "dhfr", "dhps", "mdr1", "mdr2", "arps10", "fd", "exo"
#'                  To specify a mutation put it in between quotation marks e.g. "Pfkelch13", or  c("Pfkelch13","dhfr","fd") to select several genetic mutation types
#' @param aggregate The administrative level at which we aggregate. Separate maps are created for each administrative level.
#' @param minAggregateCount The minimum count of aggregated samples. To avoid estimating on very small samples, one can set a minimum count of samples, below which the marker is not shown.
#' @param markerSize Allows adjustment of the size of markers on the map. If only one value is passed, all markers will be of that size; if two values are passed, they will be used as the min and max size of the marker, whose size will reflect the number of samples.
#' @param alleleColours A vector of R colours (specified either as R colour names, or HTML colour codes) which is used to colour the slices in the pie symbols of the plot.
#'                      All elements of the vector must be named with the name of the value (allele) they will represent.
#'                      The sequence of names will determine the order in which the values are displayed.
#'                      The name "Other" will indicate the colour used to group together samples that carry alleles other than those names. Usually, "Other" will be listed at the end of the sequence.
#'                      It no colour is specified for "Other" in the sequence, then an additional "Other"="white" element will be automatically added to the end of the sequence.
#' @param showNames If TRUE, labels are shown with the name of the aggregation unit (Province or District).
#' @param nameFontSize Allows adjustment of the Province or District label name font size (numeric), default=5.
#' @param ... Aesthetics: any of the following plot parameters,
#'            width: the width of the plot (numeric), default=15
#'            height: the height of the plot (numeric), default=15
#'            units: the units in which the width and height are expressed. Supported values: "in" (default), "cm", "mm", "px"
#'            dpi: the resolution of the plot output, expressed as dots per inch, default=300
#'            format: the file format in which the plot will be saved. Supported values: "png" (default), "pdf"
#'            legendPosition: specifies where the legend should be plotted. Supported values: "inset" (default), "separate"
#'            legendWidth: specifies how wide a fixed width legend space should be, default="NULL"
#'            legendDirection: specifies location of legend. Supported values: "vertical" (default), "horizontal"
#'            legendFontSize: specifies font size of the legend (numeric), default=4
#'            axisTitleSize: specifies axis label font size (numeric), default=1
#'
#' @export
#' @seealso Useful links:
#' \itemize{
#'  \item user guide for information on functionality of specific genetic mutation types {\url{https://genremekong.org/tools/grcmalaria-r-package-user-guide}}
#' }
#'
#' @examples \dontrun{
#' # Produce maps for all available genetic mutations for sample set called Cambodia,
#' # at province level, with a minimum of 10 samples
#' mapAlleleProportions (ctx, sampleSet="Cambodia",
#'                       mutations="ALL",
#'                       aggregate="Province", minAggregateCount=10,
#'                       showNames=TRUE, markerSize=c(10,25))
#'
#' # Produce maps of some of the most common Pfkelch13 variants (WT, C580Y, R539T and Y493H),
#' # and group together the other alleles
#' kelch13Col <- c(WT="green",C580Y="red",R539T="orange",Y493H="yellow",Other="gray")
#' mapAlleleProportions (ctx, sampleSet="Cambodia",
#'                       mutations="Pfkelch13", alleleColours=kelch13Col,
#'                       aggregate="District", minAggregateCount=10,
#'                       showNames=TRUE, markerSize=c(10,25))
#'}
#
mapAlleleProportions <- function (ctx, sampleSet, timePeriods=NULL,
                                  mutations="ALL",
                                  aggregate="Province", minAggregateCount=10,
                                  markerSize=16, alleleColours=NULL,
                                  showNames=TRUE, nameFontSize=5,
                                  ...) {

  task <- "map/alleleProp"
  params <- param.makeParameterList (ctx, task,
                                     timePeriods=timePeriods, mutations=mutations,
                                     aggregate=aggregate, minAggregateCount=minAggregateCount,
                                     markerSize=markerSize, alleleColours=alleleColours,
                                     showNames=showNames, nameFontSize=nameFontSize,
                                     ...)
  execute.executeOnSampleSet (userCtx=ctx, sampleSetName=sampleSet, task=task, params=params)
}

#############################################################
#
#' Map genetic Diversity
#'
#' This creates maps of expected heterozygosity or genetic diversity for different administrative divisions.
#' Genetic diversity can be calculated using four different measures: "maxBarcodeFreq","barcodeHet", "meanVariantHet", and "medianDistance".
#' Calculations are based on the 'genetic barcode', which consists of 101 SNPs determined for each sample. These SNPs are distributed across all nuclear chromosomes.
#' SNP positions were chosen based on their ability to differentiate populations and their power to recapitulate genetic distance.
#' If multiple measures are specified, then a map will be created for each measure.
#' For each aggregation unit, we place a marker on the map, coloured according to the genetic diversity, with a label indicating the genetic diversity.
#'
#' @param ctx The analysis context, created by intializeContext().
#' @param sampleSet The name of the sample set being used, which must have been previously created by selectSampleSet().
#' @param timePeriods Time-sequence maps can be implemented using a parameter called timePeriods parameter. When this is passed to a function, it will partition samples into time-interval plots. Time intervals parameters are available in 5 analyses: mapSampleCounts(), mapDrugResistancePrevalence(), mapMutationPrevalence(), mapAlleleProportions() and mapDiversity().
#'                    Why use time-intervals parameter:
#'                    It will produce maps with consistent geographical boundaries
#'                    Users can slice up the dataset in any specified time period
#' @param measures This can be "ALL", or any vector containing one or more of c("maxBarcodeFreq","barcodeHet", "meanVariantHet", "meanHetCallsProp", "medianDistance").
#'                 The method "maxBarcodeFreq" is a measure of loss of diversity. This measure gives the proportion of samples carrying the most common 
#'                 barcode haplotype (defined as samples with identical barcodes). The output ranges from 0-1, with a low value corresponding to a low 
#'                 proportion of samples with the most common haplotype (high diversity), and a high value corresponding to a high proportion of samples 
#'                 with the most same haplotype (low diversity).
#'                 The method "barcodeHet" is a measure of heterozygosity based on the complete barcode. This measure is the probability of two randomly 
#'                 selected samples carrying different barcodes, and is useful to detect clonal expansions.
#'                 The output ranges from 0-1, with a low value corresponding to low diversity (low probability of carrying a different barcode), 
#'                 while a high value corresponds to high diversity (high probability of a different barcode).
#'                 The method "meanVariantHet" is the average of the expected heterozygosity of the individual variants in a barcode.  For each barcode SNP, 
#'                 the expected heterozygosity is calculated using: HE=(n/(n-1))(1-sum(pi^2)), where n = the number of samples and pi = the allele frequency 
#'                 of the ith allele at this SNP. The final diversity estimate is the mean of SNP heterozygosity across all barcoding SNPs.
#'                 The output ranges from 0-1, with a low value corresponding to low diversity (low probability of a different allele), and a high value 
#'                 corresponding to high diversity (high probability of a different allele).
#'                 The method "meanHetCallsProp" is the mean of the proportion of heterozygous calls at each of the barcoding SNPs. For each barcode SNP, 
#'                 we count how many samples carry more than one allele at this SNP, and estimate the proportion, ignoring those samples that have a 
#'                 missing genotype. The final diversity estimate is the mean these proportions, across all barcoding SNPs.
#'                 The output ranges from 0-1, with a low value corresponding to low diversity, and a high value corresponding to high diversity.
#'                 The method "medianDistance" is a measure that gives the median genetic distance in geographical unit.
#'                 It is calculated using pairwise genetic distance as the proportion of SNPs differing between two samples, then taking the median 
#'                 of these in the group of samples of interest
#'                 The method "meanDistance" is a measure that gives the mean genetic distance in geographical unit, with the same method as medianDistance.
#' @param aggregate The administrative level at which we aggregate. Separate maps are created for each administrative level.
#' @param minAggregateCount The minimum count of aggregated samples. To avoid estimating on very small samples, one can set a minimum count of samples, below which the marker is not shown.
#' @param markerSize Allows adjustment of the size of markers on the map. If only one value is passed, all markers will be of that size; if two values are passed, they will be used as the min and max size of the marker, whose size will reflect the number of samples.
#' @param markerFontSize Allows adjustment of the font size shown on the markers (numeric), default=6.
#' @param markerColours The colour to indicate the level of genetic diversity. Default: "red3", other examples: "forestgreen", "cornflowerblue","darkgoldenrod3"
#' @param showNames If TRUE, labels are shown with the name of the aggregation unit (Province or District)
#' @param nameFontSize Set the font size of the geographical name labels outside the markers on the map, if showNames=TRUE, default=5
#' @param ... Aesthetics: any of the following plot parameters,
#'            width: the width of the plot (numeric), default=15
#'            height: the height of the plot (numeric), default=15
#'            units: the units in which the width and height are expressed. Supported values: "in" (default), "cm", "mm", "px"
#'            dpi: the resolution of the plot output, expressed as dots per inch, default=300
#'            format: the file format in which the plot will be saved. Supported values: "png" (default), "pdf"
#'            legendPosition: specifies where the legend should be plotted. Supported values: "inset" (default), "separate"
#'            legendWidth: specifies how wide a fixed width legend space should be, default="NULL"
#'            legendDirection: specifies location of legend. Supported values: "vertical" (default), "horizontal"
#'            legendFontSize: specifies font size of the legend (numeric), default=4
#'            axisTitleSize: specifies axis label font size (numeric), default=1
#'
#' @export
#'
#' @examples \dontrun{
#' #' #Given lists of time periods of interest, in this case calendar years
#' periods <- list(
#'                list(name="2018", type="year", start="1-Jan-2018"),
#'                list(name="2019", type="year", start="1-Jan-2019"),
#'                list(name="2020", type="year", start="1-Jan-2020")
#'            )
#'
#' # Produce genetic diversity maps for all available measures, for sampleSet "Laos",
#' # at District level, and for calendar years 2018-2020.
#' mapDiversity (ctx, sampleSet="Laos", timePeriods=periods,
#'               measures="ALL", aggregate="District", markerColours="red3",
#'               minAggregateCount=10, showNames=TRUE, markerSize=16,
#'               width=20, height=20)
#'}
#
mapDiversity <- function (ctx, sampleSet, timePeriods=NULL,
                          measures="ALL",
                          aggregate="Province", minAggregateCount=10,
                          markerSize=16, markerFontSize=6, markerColours="red3",
                          showNames=TRUE, nameFontSize=5,
                          ...) {

  task <- "map/diversity"
  #
  # If measures="ALL", we have to make separate calls, because the different measures use different datasets
  measures <- markerMap.resolveMeasureNames (ctx, "diversity", measures)
  for (measure in measures) {
      params <- param.makeParameterList (ctx, task,
                                     timePeriods=timePeriods, measures=measure, aggregate=aggregate, minAggregateCount=minAggregateCount,
                                     markerSize=markerSize, markerFontSize=markerFontSize, markerColours=markerColours, showNames=showNames, nameFontSize=nameFontSize,
                                     ...)
      execute.executeOnSampleSet (userCtx=ctx, sampleSetName=sampleSet, task=task, params=params)
  }
}


#############################################################
#
#' Map Connectedness between sites
#'
#' This creates maps of connections of haplotypes between different administrative divisions.
#' The administrative divisions are compared in a pairwise fashion, within this comparison output is produced that
#' gives either i) the proportion of sample pairs that have at least the assigned percentage of similarity (default: 1), or
#' ii) the proportion of sample pairs with a mean genetic distance of at least the meanDistanceLevels (default: 0.5).
#' Subsequently, maps are produced that show connections between the administrative divisions,
#' in which the thickness of the connection line corresponds to the proportion of either i), or ii).
#' Both similarity and mean genetic distance are calculated based on the 'genetic barcode', which consists of 101 SNPs determined for each sample.
#' If multiple levels are specified at minIdentity, and meanDistanceLevels, then multiple maps will be produced.
#'
#' @param ctx The analysis context, created by intializeContext().
#' @param sampleSet The name of the sample set being used, which must have been previously created by selectSampleSet().
#' @param measures default is "ALL", other options are "meanDistance", or "similarity".
#'                 "meanDistance" produces a map with pairwise comparisons between administrative divisions (e.g. sites) in which the thickness of the connection line corresponds to
#'                 the proportion of sample pairs with a mean genetic distance of at least the meanDistanceLevels
#'                "similarity" produces a map with pairwise comparisons between administrative divisions (e.g. sites) in which the thickness of the connection line corresponds to
#'                the proportion of sample pairs with at least the percentage of similarity given in minIdentity.
#' @param minIdentity value between 0 and 1, default is 1
#' @param meanDistanceLevels meanDistanceLevels value between 0 and 1, default is 0.5. Several thresholds can be specified at once e.g. c(1, 0.9, 0.8)
#' @param aggregate The administrative level at which we aggregate. Separate maps are created for each administrative level.
#' @param minAggregateCount The minimum count of aggregated samples. To avoid estimating on very small samples, one can set a minimum count of samples, below which the marker is not shown.
#' @param markerSize Allows adjustment of the size of markers on the map. If only one value is passed, all markers will be of that size; if two values are passed, they will be used as the min and max size of the marker, whose size will reflect the number of samples.
#' @param showNames If TRUE, labels are shown with the name of the aggregation unit (Province or District)
#' @param nameFontSize Set the font size of the geographical name labels outside the markers on the map, if showNames=TRUE, default=5
#' @param ... Aesthetics: any of the following plot parameters,
#'            width: the width of the plot (numeric), default=15
#'            height: the height of the plot (numeric), default=15
#'            units: the units in which the width and height are expressed. Supported values: "in" (default), "cm", "mm", "px"
#'            dpi: the resolution of the plot output, expressed as dots per inch, default=300
#'            format: the file format in which the plot will be saved. Supported values: "png" (default), "pdf"
#'            legendPosition: specifies where the legend should be plotted. Supported values: "inset" (default), "separate"
#'            legendWidth: specifies how wide a fixed width legend space should be, default="NULL"
#'            legendDirection: specifies location of legend. Supported values: "vertical" (default), "horizontal"
#'            legendFontSize: specifies font size of the legend (numeric), default=4
#'            axisTitleSize: specifies axis label font size (numeric), default=1
#'
#' @export
#'
#' @examples \dontrun{
#' # Map connectedness of haplotypes between districts for a similarity of 1.0, and 0.95,
#' # and a mean genetic distance of 0.5, 0.6, and 0.7
#' mapConnections (ctx, sampleSet="Laos",
#'                 measures="ALL", aggregate="District", minIdentity=c(1.0,0.95),
#'                 meanDistanceLevels=c(0.5,0.6,0.7),
#'                 minAggregateCount=10, showNames=TRUE)
#'}
#
mapConnections <- function (ctx, sampleSet,
                            measures="ALL",
                            minIdentity=1.0,
                            meanDistanceLevels=0.5,
                            aggregate="Province", minAggregateCount=10,
                            markerSize=6,
                            showNames=TRUE, nameFontSize=5,
                            ...) {

  task <- "map/connect"
  params <- param.makeParameterList (ctx, task,
                                     measures=measures, minIdentity=minIdentity, meanDistanceLevels=meanDistanceLevels,
                                     aggregate=aggregate, minAggregateCount=minAggregateCount, markerSize=markerSize,
                                     showNames=showNames, nameFontSize=nameFontSize,
                                     ...)
  execute.executeOnSampleSet (userCtx=ctx, sampleSetName=sampleSet, task=task, params=params)
}

#############################################################
#
#' Map barcode group frequencies
#'
#' Creates a map showing the frequency of groups of Plasmodium samples with identical a 'genetic barcode' within administrative divisions.
#' The samples are grouped according to identical barcode within each administrative division.
#' Barcode groups are not compared between administrative divisions. To compare barcode groups between administrative divisions please use mapClusterSharing().
#' The presence of a few large segments suggests that the location is dominated by several large groups with identical genetic barcodes.
#' Observed "black" colors in the bar or pie chart indicate the presence of many overlapping segment edges. This suggests the existence of numerous groups with low frequency. Therefore, the color black can serve as an indicator of diversity.
#'
#' @param ctx The analysis context, created by intializeContext().
#' @param sampleSet The name of the sample set being used, which must have been previously created by selectSampleSet().
#' @param type Frequency of barcode groups can be visualized either as a "bar" and/or a "pie".
#' @param aggregate The administrative level at which we aggregate. Separate maps are created for each administrative level.
#' @param minAggregateCount The minimum count of aggregated samples. To avoid estimating on very small samples, one can set a minimum count of samples, below which the marker is not shown.
#' @param markerScale Allows adjustment of the size of markers on the map, default: 0.8.
#' @param showNames If TRUE, labels are shown with the name of the aggregation unit (Province or District)
#' @param nameFontSize Set the font size of the geographical name labels outside the markers on the map, if showNames=TRUE, default=5
#' @param ... Aesthetics: any of the following plot parameters,
#'            width: the width of the plot (numeric), default=15
#'            height: the height of the plot (numeric), default=15
#'            units: the units in which the width and height are expressed. Supported values: "in" (default), "cm", "mm", "px"
#'            dpi: the resolution of the plot output, expressed as dots per inch, default=300
#'            format: the file format in which the plot will be saved. Supported values: "png" (default), "pdf"
#'            legendPosition: specifies where the legend should be plotted. Supported values: "inset" (default), "separate"
#'            legendWidth: specifies how wide a fixed width legend space should be, default="NULL"
#'            legendDirection: specifies location of legend. Supported values: "vertical" (default), "horizontal"
#'            legendFontSize: specifies font size of the legend (numeric), default=4
#'            axisTitleSize: specifies axis label font size (numeric), default=1
#'
#' @export
#'
#' @examples \dontrun{
#' # Map barcode group frequencies in both bar and pie form, for sampleSet "Laos",
#' # at both Province and District level
#' mapBarcodeFrequencies (ctx, sampleSet="Laos",
#'                        type=c("bar","pie"),
#'                        aggregate=c("Province","District"),
#'                        minAggregateCount=10, showNames=TRUE, markerScale=0.8)
#'}
#
mapBarcodeFrequencies <- function (ctx, sampleSet,
                                   type=c("bar","pie"),
                                   aggregate="Province", minAggregateCount=10,
                                   markerScale=0.8,
                                   showNames=TRUE, nameFontSize=5,
                                   ...) {

  task <- "map/barcodeFrequency"
  params <- param.makeParameterList (ctx, task,
                                     type=type, aggregate=aggregate, minAggregateCount=minAggregateCount, markerScale=markerScale,
                                     showNames=showNames, nameFontSize=nameFontSize,
                                     ...)
  execute.executeOnSampleSet (userCtx=ctx, sampleSetName=sampleSet, task=task, params=params)
}

#############################################################
#
#' Clustering
#'
#' In order to plot maps of cluster prevalence and cluster sharing, first clusters of similar genetic background need to be assigned.
#' This function partitions the Plasmodium samples into clusters with a similar genetic background, based on the genetic barcode similarity.
#' This function was built on the igraph R package (Csárdi G, 2023), which computes a graph connecting sample pairs with S greater than a minimum threshold S(min) , subsequently partitioning the graph into clusters.
#' Functions providing cluster analyses include: mapClusterSharing(), and mapClusterPrevalence().
#'
#' @param ctx The analysis context, created by intializeContext().
#' @param sampleSet The name of the sample set being used, which must have been previously created by selectSampleSet().
#' @param clusterSet Give the cluster set a name. In the example the clustering set is called "GMS".
#' @param minIdentity The minimal similarity level set for a pair of samples to be in a cluster. For example, "0.95" corresponds to at least 95 percent genetic barcode similarity.
#'                    Multiple similarity levels can be set at once, as in the example above, by putting number inside c() separated by , the default is 1.
#' @param impute To use imputed or filtered data. The default is TRUE.
#' @param clusteringMethod The clustering method. Methods available: "allNeighbours", "leiden" and "louvain". The default is "louvain".
#'                        The "allNeighbours" method clusters samples together that are above the set "minIdentity" threshold.
#'                        This method is less informative at low similarity levels, because each sample will be assigned to a single cluster.
#'                        The "louvain" and "leiden" methods are preferred. They use an algorithm to identify clusters within a network that 
#'                        are strongly connected to each other, and more weakly connected to other clusters.
#' @param minClusterSize To avoid creating very small clusters, one can set a minimum cluster size. The default is 10 samples.
#' @param resolution A tuning parameter for the "leiden" and "louvain" methods. 
#'                        Values are between 0 and 1; lowering the value produces larger, "looser" clusters. Default is 1. 
#' @param objectiveFunction A parameter for the "leiden" method only. Possible values: "CPM" (Constant Potts Model) or "modularity". Default is "CPM".
#' @param beta A tuning parameter for the "leiden" method only. Affects the randomness in the algorithm. Default is 0.01.
#'
#' @export
#'
#' @examples \dontrun{
#' ## Find clusters of similar genetic background ##
#'    findClusters(ctx, sampleSet="EBKK", clusterSet = "GMS",
#'                 minIdentity = c(0.95, 0.80), impute = TRUE,
#'                 clusteringMethod = "louvain")
#'}
#
findClusters <- function (ctx, sampleSet, clusterSet,
                          minIdentity=1.0, impute=TRUE,
                          clusteringMethod="louvain", minClusterSize=10, 
                          resolution = 1, objectiveFunction = "CPM", beta = 0.01) {
  # Construct a list of arguments so we can use the params.getArgParameter utility function
  args <- list(clusterSet=clusterSet, minIdentity=minIdentity, impute=impute, clusteringMethod=clusteringMethod, minClusterSize=minClusterSize)
  p <- new.env ()
  p$cluster.clusterSet.name <- param.getArgParameter (args, "clusterSet")
  p$cluster.identity.min    <- param.getArgParameter (args, "minIdentity",      type="numeric", multiValue=TRUE, defaultValue=1.0)
  p$cluster.impute          <- param.getArgParameter (args, "impute",           type="logical", defaultValue=TRUE)
  p$cluster.minSize         <- param.getArgParameter (args, "minClusterSize",   type="integer", defaultValue=10)
  p$cluster.method          <- param.getArgParameter (args, "clusteringMethod", defaultValue="louvain", validValues=c("louvain","leiden","allNeighbours"))
  p$cluster.resolution      <- param.getArgParameter (args, "resolution",       type="numeric", defaultValue=1.0)
  p$cluster.objFunction     <- param.getArgParameter (args, "objectiveFunction", defaultValue="CPM", validValues=c("CPM","modularity"))
  p$cluster.beta            <- param.getArgParameter (args, "beta",             type="numeric", defaultValue=0.01)
  
  cluster.findClusters (ctx, sampleSetName=sampleSet, params=p)
}

#############################################################
#
#' Plot graphs of barcode identity, using cluster
#'
#' Creates graphs of cluster networks within a sampleset of Plasmodium samples.
#' The graphs are based on a cluster analysis, which needs to be produced with the findClusters function in advance.
#' The findClusters function partitions the Plasmodium samples into clusters with a similar genetic background, based on the 'genetic barcode' similarity.
#'
#' @param ctx The analysis context, created by intializeContext()
#' @param sampleSet The name of the sample set being used, which must have been previously created by selectSampleSet()
#' @param clusterSet The name of the clustering set, defined in findClusters()
#' @param graphLayout TBD, the default="fr"
#' @param weightPower TBD, the default is 2
#' @param ... Aesthetics: any of the following plot parameters,
#'            width: the width of the plot (numeric), default=15
#'            height: the height of the plot (numeric), default=15
#'            units: the units in which the width and height are expressed. Supported values: "in" (default), "cm", "mm", "px"
#'            dpi: the resolution of the plot output, expressed as dots per inch, default=300
#'            format: the file format in which the plot will be saved. Supported values: "png" (default), "pdf"
#'            legendPosition: specifies where the legend should be plotted. Supported values: "inset" (default), "separate"
#'            legendWidth: specifies how wide a fixed width legend space should be, default="NULL"
#'            legendDirection: specifies location of legend. Supported values: "vertical" (default), "horizontal"
#'            legendFontSize: specifies font size of the legend (numeric), default=4
#'            axisTitleSize: specifies axis label font size (numeric), default=1
#'
#' @return Produces network graphs of clusters assigned in findClusters function.
#' @export
#'
#' @examples \dontrun{
#' #Plot cluster network graph for sampleset "Laos", and clusterset "LAclust"
#' plotClusterGraph (ctx, sampleSet="Laos", clusterSet="LAclust",
#'                   graphLayout="fr", weightPower=2)
#'}
#
plotClusterGraph <- function (ctx, sampleSet, clusterSet,
                              graphLayout="fr", weightPower=2,
                              ...) {

  task <- "graph"
  params <- param.makeParameterList (ctx, task,
                                     clusterSet=clusterSet, graphLayout=graphLayout, weightPower=weightPower,
                                     ...)
  execute.executeOnSampleSet (userCtx=ctx, sampleSetName=sampleSet, task=task, params=params)
}

#############################################################
#
#' Cluster Group Sharing Maps
#'
#' This creates a map showing the proportions of assigned clusters for different administrative divisions.
#' The different colours represent different haplotype groups.
#' White represents samples that do not belong to a cluster because they did not meet the specified parameter criteria.
#' The maps are based on the cluster analysis context, which needs to be produced with the findClusters function in advance.
#'
#' @param ctx The analysis context, created by intializeContext().
#' @param sampleSet The name of the sample set being used, which must have been previously created by selectSampleSet().
#' @param clusterSet The name of the cluster set, defined in findClusters().
#' @param type Frequency of clusters can be visualized either as a "bar" and/or a "pie". The default is c("bar", "pie").
#' @param aggregate The administrative level at which we aggregate. Separate maps are created for each administrative level.
#' @param minAggregateCount The minimum count of aggregated samples. To avoid estimating on very small samples, one can set a minimum count of samples, below which the marker is not shown.
#' @param markerScale Allows adjustment of the size of markers on the map, default: 0.8.
#' @param showNames If TRUE, labels are shown with the name of the aggregation unit (Province or District)
#' @param nameFontSize Set the font size of the geographical name labels outside the markers on the map, if showNames=TRUE, default=5
#' @param ... Aesthetics: any of the following plot parameters,
#'            width: the width of the plot (numeric), default=15
#'            height: the height of the plot (numeric), default=15
#'            units: the units in which the width and height are expressed. Supported values: "in" (default), "cm", "mm", "px"
#'            dpi: the resolution of the plot output, expressed as dots per inch, default=300
#'            format: the file format in which the plot will be saved. Supported values: "png" (default), "pdf"
#'            legendPosition: specifies where the legend should be plotted. Supported values: "inset" (default), "separate"
#'            legendWidth: specifies how wide a fixed width legend space should be, default="NULL"
#'            legendDirection: specifies location of legend. Supported values: "vertical" (default), "horizontal"
#'            legendFontSize: specifies font size of the legend (numeric), default=4
#'            axisTitleSize: specifies axis label font size (numeric), default=1
#'
#' @return Maps with cluster frequencies across administrative divisions
#' @export
#'
#' @examples \dontrun{
#' # Map cluster frequencies in both bar and pie form, for sampleset "Laos",
#' # and clusterset "LAclust" at both Province and District level,
#' mapClusterSharing (ctx, sampleSet="Laos", clusterSet = "LAclust",
#'                    type=c("bar", "pie"),
#'                    aggregate=c("Province","District"),
#'                    minAggregateCount=5, showNames=TRUE, markerScale=0.8)
#'}
#
mapClusterSharing <- function (ctx, sampleSet, clusterSet,
                               type=c("bar","pie"),
                               aggregate="Province", minAggregateCount=5,
                               markerScale=0.8,
                               showNames=TRUE, nameFontSize=5,
                               ...) {

  task <- "map/clusterSharing"
  params <- param.makeParameterList (ctx, task,
                                     clusterSet=clusterSet, type=type, aggregate=aggregate, minAggregateCount=minAggregateCount, markerScale=markerScale,
                                     showNames=showNames, nameFontSize=nameFontSize,
                                     ...)
  execute.executeOnSampleSet (userCtx=ctx, sampleSetName=sampleSet, task=task, params=params)
}

#############################################################
#
#' Cluster Prevalence Maps
#'
#' This function produces network maps showing the prevalence of clusters.
#' For each assigned cluster a map is created, showing the prevalence and overlap in prevalence of that cluster across administrative divisions.
#' Each map also contains information on sample count, drug-resistance, and mutation frequencies.
#' The maps are based on the cluster analysis, which needs to be produced with the findClusters function in advance.
#' The prevalence of the assessed cluster is displayed in the nodes.
#' The overlap in prevalence is displayed in the edges (i.e. lines).
#' The definition of prevalence overlap is the proportion of sample pairs that have at least the assigned "minIdentity" of genetic barcode similarity in the findCluster function.
#'
#' @param ctx The analysis context, created by intializeContext().
#' @param sampleSet The name of the sample set being used, which must have been previously created by selectSampleSet().
#' @param clusterSet The name of the cluster set, defined in findClusters().
#' @param aggregate The administrative level at which we aggregate. Separate maps are created for each administrative level.
#' @param minAggregateCount The minimum count of aggregated samples. To avoid estimating on very small samples, one can set a minimum count of samples, below which the marker is not shown.
#' @param markerSize Allows adjustment of the size of markers on the map. If only one value is passed, all markers will be of that size; if two values are passed, they will be used as the min and max size of the marker, whose size will reflect the number of samples.
#' @param markerFontSize Allows adjustment of the font size shown on the markers (numeric), default=6.
#' @param showNames If TRUE, labels are shown with the name of the aggregation unit (Province or District)
#' @param nameFontSize Set the font size of the geographical name labels outside the markers on the map, if showNames=TRUE, default=5
#' @param ... Aesthetics: any of the following plot parameters,
#'            width: the width of the plot (numeric), default=15
#'            height: the height of the plot (numeric), default=15
#'            units: the units in which the width and height are expressed. Supported values: "in" (default), "cm", "mm", "px"
#'            dpi: the resolution of the plot output, expressed as dots per inch, default=300
#'            format: the file format in which the plot will be saved. Supported values: "png" (default), "pdf"
#'            legendPosition: specifies where the legend should be plotted. Supported values: "inset" (default), "separate"
#'            legendWidth: specifies how wide a fixed width legend space should be, default="NULL"
#'            legendDirection: specifies location of legend. Supported values: "vertical" (default), "horizontal"
#'            legendFontSize: specifies font size of the legend (numeric), default=4
#'            axisTitleSize: specifies axis label font size (numeric), default=1
#'
#' @return For each assigned cluster a map with cluster prevalence and overlap in cluster prevalence across administrative divisions
#' @export
#'
#' @examples \dontrun{
#' # Produce a network and prevalence map for each assigned cluster,
#' # for sampleSet "Laos", and clusterSet "LAclust" at both a Province and District level,
#' mapClusterPrevalence (ctx, sampleSet="Laos", clusterSet = "LAclust",
#'                       aggregate=c("Province","District"),
#'                       minAggregateCount=5, showNames=TRUE)
#'}
#
mapClusterPrevalence <- function (ctx, sampleSet, clusterSet,
                                  aggregate="Province", minAggregateCount=5,
                                  markerSize=16, markerFontSize=6,
                                  showNames=TRUE, nameFontSize=5,
                                  ...) {

  task <- "map/clusterPrevalence"
  params <- param.makeParameterList (ctx, task,
                                     clusterSet=clusterSet, aggregate=aggregate, minAggregateCount=minAggregateCount, markerSize=markerSize,
                                     markerFontSize=markerFontSize, showNames=showNames, nameFontSize=nameFontSize,
                                     ...)
  execute.executeOnSampleSet (userCtx=ctx, sampleSetName=sampleSet, task=task, params=params)
}

#############################################################
#
#' Set to a user-defined colour palette
#'
#' Users can define their own colour palette instead of using the default colours.
#' This function needs to be executed before sample selection, and after running initializeContext(), otherwise the defined colour palette will not be applied.
#'
#' @param ctx The analysis context, created by intializeContext().
#' @param palette The desired palette, specified as separate HTML colour codes
#'
#' @export
#'
#' @examples \dontrun{
#' # To view default 25-colour palette
#' ctx$config$defaultPalette
#'
#' # Set colours
#' setColourPalette(ctx, palette=c("red","lightblue","#93c4d2"))
#'
#' # then select sample set
#' selectSampleSet(ctx, sampleSet="SamplesetName",
#'                 select=list(
#'                            list(field="Species", values=c("Pf"))
#'                 ))
#'
#' # then proceed with the analysis
#'}
#
setColourPalette <- function (ctx, palette) {
  graphics.setColourPalette (ctx, palette)
}

#############################################################
#
#' Reset to the default colour palette, removing user-defined palettes
#'
#' Return to the default 25-colour palette of grcMalaria.
#'
#' @param ctx The analysis context, created by intializeContext().
#'
#' @return Return to the default colour palette of grcMalaria.
#' @export
#'
#' @examples \dontrun{
#' # return to the default 25-colour palette:
#' resetColourPalette(ctx)
#'
#' # To view default 25-colour palette
#' ctx$config$defaultPalette
#'}
#
resetColourPalette <- function (ctx) {
  graphics.resetColourPalette (ctx)
}

#############################################################
#
#' Graphical Attributes for Principal Component Analysis (PCA) and Neighbor-Joining (NJ)-tree
#'
#' This function allows for specification of the graphical display of samples based on their metadata,
#' that can subsequently be used to plot the population structure with a principal component analysis (plotPrincipalComponents()),
#' or creation of a neighbour-joining tree (plotTree()).
#' Specification based on metadata (i.e. attributes) is done by loading an Excel spreadsheet into R.
#' The spreadsheet and instructions can be found in the grcMalaria user-guide, see link in 'See also' section.
#'
#' @param ctx The analysis context, created by intializeContext()
#' @param name The name of the graphical attribute, this can be any metadata category of interest
#'             Once defined, graphical attributes can be applied by referencing their name
#'             Examples of names are: 'Countries', 'Provinces', 'Kelch13', etc
#' @param field The corresponding column-name in the original data file of the specified graphical attributes
#' @param file The file name of the Excel spreadsheet where all graphical attributes are specified
#' @param sheet The sheet name of the Excel spreadsheet where the to parameter 'name' corresponding graphical attribute is specified
#'
#' @return The analysis context, augmented with the new sample set
#' @export
#'
#' @seealso grcMalaria user-guide:
#' \itemize{
#'  \item {\url{https://genremekong.org/tools/grcmalaria-r-package-user-guide}}
#' }
#'
#' @examples \dontrun{
#' #Load Excel spreadsheet into R with specified graphical attributes
#' gaFile <- "/path/to/my/file/GenRe-GraphicAttributes.xlsx"
#'
#' #Load graphical attributes for province, kelch13, and plasmepsin2/3 amplification status
#' ctx <- loadGraphicAttributes (ctx, name="province", field="AdmDiv1",
#'                               file=gaFile, sheet="Provinces")
#' ctx <- loadGraphicAttributes (ctx, name="k13", field="Kelch",
#'                               file=gaFile, sheet="Kelch13")
#' ctx <- loadGraphicAttributes (ctx, name="pm23", field="Plasmepsin2/3",
#'                               file=gaFile, sheet="Plasmepsin23")
#' ctx <- loadGraphicAttributes (ctx, name="pm23-noColour", field="Plasmepsin2/3",
#'                               file=gaFile, sheet="Plasmepsin23-noColour")
#'}
#
loadGraphicAttributes <- function (ctx, name, field, file, sheet) {
  graphics.loadAttributes (ctx, name, field, file, sheet)
}

#############################################################
#
#' Principal Component Analyses (PCA)
#'
#' Principal Component Analysis (PCA) or Principal Coordinate Analysis (PCoA) are methods of reducing complexity of the data (in this case distances of 101-SNP genetic barcodes) down to two dimensions, while preserving as much information contained in the original data as possible.
#' This was done by creating principal components that explain most of the observed variance in the dataset.
#' PCoA is the default method in the grcMalaria but other PCA methods are available.
#' One plot is created for each specified graphical attribute in the parameter 'plots', in which the samples are coloured according to user input.
#' Graphical display of attributes needs to be specified before using this function.
#' The PC (Principal Components) are arranged in descending order of explained variance, with PC1 axis showing the largest variation and PC2 axis showing the second the largest variation and so on.
#' Samples that are close to each other means they are similar to each other. Well separated clusters means there are differences between samples.
#'
#' @param ctx The analysis context, created by intializeContext()
#' @param sampleSet The name of the sample set being used, which must have been previously created by selectSampleSet()
#' @param type The type of method used for the principal component analysis
#'             "PCoA", performs multidimensional scaling on a distance matrix, in which each sample is a variable
#'             "nipals", "bpca" are methods applied to the set of genetic barcode genotypes (each barcode position is a variable).
#'             The default is "PCoA"
#' @param plots The list of attributes for which plots will be created.
#' @param impute To use imputed or filtered data. The default is TRUE.
#' @param ... Aesthetics: any of the following plot parameters,
#'            width: the width of the plot (numeric), default=15
#'            height: the height of the plot (numeric), default=15
#'            units: the units in which the width and height are expressed. Supported values: "in" (default), "cm", "mm", "px"
#'            dpi: the resolution of the plot output, expressed as dots per inch, default=300
#'            format: the file format in which the plot will be saved. Supported values: "png" (default), "pdf"
#'            legendPosition: specifies where the legend should be plotted. Supported values: "inset" (default), "separate"
#'            legendWidth: specifies how wide a fixed width legend space should be, default="NULL"
#'            legendDirection: specifies location of legend. Supported values: "vertical" (default), "horizontal"
#'            legendFontSize: specifies font size of the legend (numeric), default=4
#'            axisTitleSize: specifies axis label font size (numeric), default=1
#'
#' @return Principal component plots and .tab data files.
#' @export
#'
#' @examples \dontrun{
#' # Perform a principal component analysis for province, kelch13, plasmepsin2/3 amplification
#' # status, and kelch13 in combination with plasmepsin2/3 amplification status,
#' # using the "PCoA" method on sampleSet "Laos"
#' plotPrincipalComponents (ctx, sampleSet="Laos", type="PCoA",
#'     plots=list(
#'               list(name="ByProvince",       attributes="province"),
#'               list(name="ByKelch13",        attributes="k13"),
#'               list(name="ByPm23",           attributes="pm23"),
#'               list(name="ByKelch13AndPm23", attributes=c("k13","pm23-noColour"))
#'            ),
#'     width=15, height=10)
#'}
#
plotPrincipalComponents <- function (ctx, sampleSet,
                                     type="PCoA", plots, impute=FALSE,
                                     ...) {

  # Construct a list of arguments so we can use the params.getArgParameter utility function
  args <- list(type=type, plots=plots)
  pcaType  <- param.getArgParameter (args, "type", defaultValue="PCoA", validValues=c("PCoA", "nipals", "bpca"))
  task <- paste("pca", pcaType, sep="/")
  params <- param.makeParameterList (ctx, task,
                                     plots=plots,
                                     impute=impute,
                                     ...)

  execute.executeOnSampleSet (userCtx=ctx, sampleSetName=sampleSet, task=task, params=params)
}

#############################################################
#
#' Plot Trees
#'
#' Neighbor-joining (NJ) tree is a another method to show how samples are related based on the distances between barcodes.
#' It is not a phylogenetic tree. The NJ tree is unrooted, which means: 1) It clusters together related sequences using the distance matrix,
#' the tree has no orientation, there is no assumption about ancestry.
#' For each specified attribute, a plot is created for the phylogenetic tree, in which the nodes (i.e. samples) are coloured according to that attribute.
#' Graphical display of attributes needs to be specified before using this function.
#' In addition to the plot, a newick file is generated, which can be used to upload in tree visualisation tools (for example https://microreact.org/)
#'
#' @param ctx The analysis context, created by intializeContext().
#' @param sampleSet The name of the sample set being used, which must have been previously created by selectSampleSet().
#' @param type The method used to create the tree.Currently the only option "njt" which creates a neighbour-joining tree.
#' @param plots The list of attributes for which plots will be created. See example below.
#' @param impute To use imputed or filtered data. The default is TRUE.
#' @param ... Aesthetics: any of the following plot parameters,
#'            width: the width of the plot (numeric), default=15
#'            height: the height of the plot (numeric), default=15
#'            units: the units in which the width and height are expressed. Supported values: "in" (default), "cm", "mm", "px"
#'            dpi: the resolution of the plot output, expressed as dots per inch, default=300
#'            format: the file format in which the plot will be saved. Supported values: "png" (default), "pdf"
#'            legendPosition: specifies where the legend should be plotted. Supported values: "inset" (default), "separate"
#'            legendWidth: specifies how wide a fixed width legend space should be, default="NULL"
#'            legendDirection: specifies location of legend. Supported values: "vertical" (default), "horizontal"
#'            legendFontSize: specifies font size of the legend (numeric), default=4
#'            axisTitleSize: specifies axis label font size (numeric), default=1
#'
#' @return Phylogenetic tree plots, .tab data files, and a .newick tree file.
#' @export
#'
#' @examples \dontrun{
#' # Produce a neighbour-joining tree for sampleSet "Laos"
#' # with plots in which the samples (i.e. nodes) are coloured according to i) province, ii) kelch13,
#' # iii) plasmepsin2/3 amplification status, iv) kelch13 and plasmepsin2/3 amplification status
#' plotTree (ctx, sampleSet="Laos", type="njt",
#'           plots=list(
#'                     list(name="ByProvince",       attributes="province"),
#'                     list(name="ByKelch13",        attributes="k13"),
#'                     list(name="ByPm23",           attributes="pm23"),
#'                     list(name="ByKelch13AndPm23", attributes=c("k13","pm23-noColour"))
#'           ),
#'           width=15, height=10)
#'}
#
plotTree <- function (ctx, sampleSet,
                      type="njt", plots, impute=FALSE,
                      ...) {

  # Construct a list of arguments so we can use the params.getArgParameter utility function
  args <- list(type=type, plots=plots)
  treeType  <- param.getArgParameter (args, "type", defaultValue="njt", validValues=c("njt"))
  task <- paste("tree", treeType, sep="/")
  params <- param.makeParameterList (ctx, task,
                                     plots=plots,
                                     impute=impute,
                                     ...)
  execute.executeOnSampleSet (userCtx=ctx, sampleSetName=sampleSet, task=task, params=params)
}
#
#
#
#############################################################
