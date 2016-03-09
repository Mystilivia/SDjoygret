#' First function
#'
#' This function is just a test for R package creation.
#' @param x Must be numeric
#' @keywords Test
#' @export
#' @examples
#' testing_function_01(1)

testing_function_01 <- function(x){
  if (is.numeric(x)==FALSE){
    stop("Please give me numeric to eat !")
  }
  if(x==1){
    print("Yeah ! you typed 1")
    result.step1 <- x*5
    result.step2 <- x/2
    list.result <- list(result.step1, result.step2)
    return(list.result)
  } else {
    print("I don't get it...")
  }
}

#' xcms_orbi_GRT
#'
#' This function simplify workflow from peak grouping and retention time correction
#' to PCA analysis. Function will perform group -> retcor -> group -> retcor -> group
#' and return EIC for standards across all samples and extract standards from peaks table,
#' giving also ppm deviation for each std for diagnostic.
#' @param File_list File list to pass to xcmsSet method = 'centWave', ppm=7, peakwidth=c(4,20), snthresh=10, prefilter=c(4,1000), mzdiff=-0.001, fitgauss=F, nSlaves = 4
#' @param Results.dir.name Name of the subfolder to store results
#' @param bw_param Vector for bw settings to use in 1, 2 and 3 iteration : c(1, 2, 3)
#' @param mzwid_param mzwid parameter to use.
#' @param minfrac_param minfrac parameter to use
#' @param profStep_param profStep parameter to use
#' @param Sample.Metadata Dataframe with samples names and metadata
#' @param Grouping.factor Col.number of Sample.Metadata to use for ACP colors (can be a vector)
#' @keywords xcms, orbitrap
#' @usage xcms_orbi_GRT(xcms_set_obj, Results.dir.name="Default", bw_param=c(25, 10, 0.7), mzwid_param=0.005, minfrac_param=0.25, profStep_param=0.8)
#' xcms_orbi_GRT()

xcms_orbi_GRT <- function(File_list,
                          Results.dir.name="Default",
                          bw_param=c(15, 8, 0.8),
                          mzwid_param=0.005,
                          minfrac_param=0.7,
                          profStep_param=0.5,
                          Sample.Metadata,
                          Grouping.factor=1){
  ## Package requirement
  require("xcms")
  require("ropls")

  ## Create directory and path
  Results.path.root <- paste0("./",Results.dir.name,"/")
  Results.path.rtgraph <- paste0(Results.path.root, "/RetCor_Dev/")
  Results.path.pca <- paste0(Results.path.root, "/PCA/")
  dir.create(Results.path.root, showWarnings = F)
  dir.create(Results.path.rtgraph, showWarnings = F)
  dir.create(Results.path.pca, showWarnings = F)

  xset.default <- xcmsSet(File_list, method = 'centWave', ppm=7, peakwidth=c(4,20), snthresh=10, prefilter=c(4,1000), mzdiff=-0.001, fitgauss=F, nSlaves = 4)
  xset.group <- xcms::group(xset.default, method="density", bw=bw_param[1], mzwid=mzwid_param, minfrac=mzwid_param)
  xset.2 <- xcms::retcor(xset.group, method="obiwarp", profStep=profStep_param, plottype="deviation")
  dev.copy(png, paste0(Results.path.rtgraph, "RetCor_01.png"), h=800, w=1600)
  dev.off()
  xset.group2 <- xcms::group(xset.2, method="density", bw=bw_param[2], mzwid=mzwid_param, minfrac=mzwid_param)
  xset.3 <- xcms::retcor(xset.group2, method="obiwarp", profStep=profStep_param, plottype="deviation")
  dev.copy(png, paste0(Results.path.rtgraph, "RetCor_02.png"), h=800, w=1600)
  dev.off()
  xset.group.4 <- xcms::group(xset.3, method="density", bw=bw_param[3], mzwid=mzwid_param, minfrac=mzwid_param)
  xset.filled <- xcms::fillPeaks(xset.group.4)
  print(xset.filled) ## Output results

  ## Generate Data.matrix, Samples and Variables metadata
  write.table(peakTable(xset.filled), file=paste0(Results.path.root, "Peak_Table.csv"), sep=";", col.names=NA)
  Sample.Metadata.D <- data.frame(row.names=Sample.Metadata[,1], Sample.Metadata[2:ncol(Sample.Metadata)])
  Data <- list()
  Data[[1]] <- t(peakTable(xset.filled)[(ncol(peakTable(xset.filled))-nrow(xset.filled@phenoData)+1):ncol(peakTable(xset.filled))])
  names(Data[1]) <- "Datamatrix"
  Data[[2]] <- subset(Sample.Metadata.D, rownames(Sample.Metadata.D) %in% rownames(xset.filled@phenoData))
  names(Data[2]) <- "Sample.metadata"
  Data[[3]] <- peakTable(xset.filled)[1:(ncol(peakTable(xset.filled))-nrow(xset.filled@phenoData))]
  names(Data[3]) <- "Variable.metadata"

  print(summary(Data)) ## Debug
  ## internal STD 1 133.1062 rt 60  Leucine-5,5,5,-d3
  ## internal STD 2 206.1014 rt 151 L-tryptophan-2,3,3-d3
  ## internal STD 3 179.0874 rt 375 Indole-2,4,5,6,7-d5-3-acetic acid
  ## internal STD 4 281.3265 rt 566 1,14-tetradecanedioic-d24 acid
  group.id.p <- as.numeric(c(rownames(subset(Data[[3]], Data[[3]][, "mz"] >= min(xcms:::ppmDev(133.1062, 5)) & Data[[3]][, "mz"] <= max(xcms:::ppmDev(133.1062, 5)))),
                             rownames(subset(Data[[3]], Data[[3]][, "mz"] >= min(xcms:::ppmDev(206.1014, 5)) & Data[[3]][, "mz"] <= max(xcms:::ppmDev(206.1014, 5)))),
                             rownames(subset(Data[[3]], Data[[3]][, "mz"] >= min(xcms:::ppmDev(179.0874, 5)) & Data[[3]][, "mz"] <= max(xcms:::ppmDev(179.0874, 5)))),
                             rownames(subset(Data[[3]], Data[[3]][, "mz"] >= min(xcms:::ppmDev(281.3265, 5)) & Data[[3]][, "mz"] <= max(xcms:::ppmDev(281.3265, 5))))))

  temp.eic.r <- getEIC(xset.filled, groupidx=group.id.p, rt="raw")
  temp.eic.c <- getEIC(xset.filled, groupidx=group.id.p, rt="corrected")
  temp.n <- length(group.id.p)
  par(mfrow=c(temp.n,2))
  for (i in 1:length(group.id.p)){
    plot(temp.eic.r, xset.filled, groupidx=i, main="RAW")
    plot(temp.eic.c, xset.filled, groupidx=i, main="Corrected")
  }
  dev.copy(png, paste0(Results.path.root, "EIC.STD.peaks.png"), h=1400, w=1000)
  dev.off()
  par(mfrow=c(1,1))

  ## Write table with STDs only and ppm deviation
  STDs.results <- subset(Data[[3]], rownames(Data[[3]]) %in% group.id.p)
  STDs.results$ppm <- (STDs.results$mzmax - STDs.results$mzmin)/(STDs.results$mz/1000000)
  STDs.results
  write.table(STDs.results, file=paste0(Results.path.root, "Results_STDs.csv"), sep=";", col.names = NA)

  ## Perform PCA
  library(ropls)
  ACP.results.list <- list()
  ACP.results.list[["ACP"]] <- opls(Data[[1]], predI=NA, plotL=F)
  for (i in 1:length(ACP.results.list)){
    Data.pca <- ACP.results.list[[i]]
    png(filename=paste0(Results.path.pca, names(ACP.results.list[i]),".png"), width=800, height=1200, units="px", res=150)
    par(mfrow=c(3,2))
    plot(Data.pca, typeVc="overview", parDevNewL=F)
    plot(Data.pca, typeVc="x-loading", parDevNewL=F)
    plot(Data.pca, typeVc="x-score", parDevNewL=F)
    plot(Data.pca, typeVc="outlier", parDevNewL=F)
    plot(Data.pca, typeVc="correlation", parDevNewL=F)
    dev.off()
    rm(Data.pca)
  }
  x <- -4
  y <- -7
  png(filename=paste0(Results.path.pca, "ACP_Ellipses.png"), width=900, height=900, units="px", res=100)
  par(mfrow=c(2,2))
  for (i in Grouping.factor){
    Data.pca <- ACP.results.list[[1]]
    temp.factor <- Data[[2]][,i]
    temp.factor.names <- names(Data[[2]][i])
    plot(Data.pca, typeVc="x-score", parAsColFcVn=addNA(as.factor(temp.factor)), parEllipses=F, parDevNewL=F)
    text(x,y,temp.factor.names)
  }
  dev.off()
  return(xset.filled) ## Output results
}


#' xcms_orbi_A
#'
#' This function perform the first steps of metabolomic analysis : xcmsSet and
#' two iteration with group and retcor to end with fillpeaks. The resulting output
#' is a xcms set object.
#' @param File_list File list to pass to xcmsSet method = 'centWave', ppm=7, peakwidth=c(4,20), snthresh=10, prefilter=c(4,1000), mzdiff=-0.001, fitgauss=F, nSlaves = 4
#' @param xcmsSet_param xcmsSet parameters for centwave method : c(ppm, min_peakwidth, max_peakwidth, snthresh, prefilter_scan, prefilter_value, mzdiff, fitgauss, nSlaves)
#' @param Results.dir.name Name of the subfolder to store results
#' @param bw_param Vector for bw settings to use in 1, 2 and 3 iteration : c(1, 2, 3)
#' @param mzwid_param mzwid parameter to use.
#' @param minfrac_param minfrac parameter to use
#' @param profStep_param profStep parameter to use
#' @keywords xcms, orbitrap
#' @usage xcms_orbi_GRT(xcms_set_obj, Results.dir.name="Default", bw_param=c(25, 10, 0.7), mzwid_param=0.005, minfrac_param=0.25, profStep_param=0.8)
#' xcms_orbi_GRT()
#' @export

xcms_orbi_A <- function(File_list,
                        xcmsSet_param=c(7, 4, 20, 10, 4, 1000, -0.001, "F", 4),
                          Results.dir.name="Default",
                          bw_param=c(15, 8, 0.8),
                          mzwid_param=0.005,
                          minfrac_param=0.7,
                          profStep_param=0.5){
  ## Package requirement
  require("xcms")

  ## Create directory and path
  Results.path.root <- paste0("./",Results.dir.name,"/")
  dir.create(Results.path.root, showWarnings = F)
  xset.default <- xcmsSet(File_list, method = 'centWave',
                          ppm=as.numeric(xcmsSet_param[1]),
                          peakwidth=as.numeric(c(xcmsSet_param[2],xcmsSet_param[3])),
                          snthresh=as.numeric(xcmsSet_param[4]),
                          prefilter=as.numeric(c(xcmsSet_param[5],xcmsSet_param[6])),
                          mzdiff=as.numeric(xcmsSet_param[7]),
                          fitgauss=xcmsSet_param[8],
                          nSlaves = as.numeric(xcmsSet_param[9]))
  xset.group <- xcms::group(xset.default, method="density", bw=bw_param[1], mzwid=mzwid_param, minfrac=mzwid_param)
  xset.2 <- xcms::retcor(xset.group, method="obiwarp", profStep=profStep_param, plottype="deviation")
  dev.copy(png, paste0(Results.path.rtgraph, "RetCor_01.png"), h=800, w=1600)
  dev.off()
  xset.group2 <- xcms::group(xset.2, method="density", bw=bw_param[2], mzwid=mzwid_param, minfrac=mzwid_param)
  xset.3 <- xcms::retcor(xset.group2, method="obiwarp", profStep=profStep_param, plottype="deviation")
  dev.copy(png, paste0(Results.path.rtgraph, "RetCor_02.png"), h=800, w=1600)
  dev.off()
  xset.group.4 <- xcms::group(xset.3, method="density", bw=bw_param[3], mzwid=mzwid_param, minfrac=mzwid_param)
  xset.filled <- xcms::fillPeaks(xset.group.4)
  return(xset.filled) ## Output results
}




