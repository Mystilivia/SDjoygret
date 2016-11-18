#' First function
#'
#' This function is just a test for R package creation.
#' @param x Must be numeric
#' @keywords Test
#' testing_function_01()
#' @export
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
                          Results.dir.name = "Default",
                          bw_param = c(15, 8, 0.8),
                          mzwid_param = 0.005,
                          minfrac_param = 0.7,
                          profStep_param = 0.5,
                          Sample.Metadata,
                          Grouping.factor = 1){
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
#' is a xcms set object. Parameters are saved in a table.
#' @param File_list        File list to pass to xcmsSet method = 'centWave', ppm=7, peakwidth=c(4,20), snthresh=10, prefilter=c(4,1000), mzdiff=0.015, fitgauss=F, nSlaves = 4, phenoData = Sample.Metadata)
#' If phenoData is a dataframe with at lest a 'class' column. Could be any number of columns but in files sample order. If 'injectionorder' and 'batch" columns are
#' presents, those will be use to order en color QCs plot.
#' @param xcmsSet_param    xcmsSet parameters (can be any parameters formated like this : c(method="centWave", "ppm=7)).
#' @param Results.dir.name Name of the subfolder to store results.
#' @param bw_param         Vector for bw settings to use in 1, 2 and 3 iteration : c(1, 2, 3).
#' @param mzwid_param      mzwid parameter to use.
#' @param minfrac_param    minfrac parameter to use.
#' @param profStep_param   profStep parameter to use.
#' @param STDs_data        A dataframe of m/z to plot with at least one "mz" named column.
#' @param STDs_EIC         Should EIC of ions list be plotted.
#' @param QCs_Graph        Logical to determine of QCs graphs needs to be saved.
#' @keywords xcms, orbitrap
#' @usage xcms_orbi_A()
#' xcms_orbi_A()
#' @export

xcms_orbi_A <- function(File_list,
                        xcmsSet_param    = list(method="centWave",
                                                ppm=7,
                                                peakwidth=c(4,20),
                                                snthresh=10,
                                                prefilter=c(4,10000),
                                                mzdiff=-0.015,
                                                fitgauss=FALSE,
                                                nSlaves=4,
                                                phenoData = Sample.Metadata),
                        Results.dir.name = Results.path,
                        bw_param         = c(15, 8, 0.8),
                        mzwid_param      = 0.015,
                        minfrac_param    = 0.7,
                        profStep_param   = 0.5,
                        STDs_data        = NULL,
                        STDs_EIC         = FALSE,
                        QCs_Graph        = FALSE) {
  ## Package requirement
  require("xcms")

  ## Create directory and path

  i <- 1
  Results.path.root <- paste0(Results.path, "XCMS_Result_", i, "/")
  while (dir.exists(Results.path.root)==TRUE) {
    i <- i+1
    Results.path.root <- paste0(Results.path, "XCMS_Result_", i, "/")
  }
  dir.create(path = Results.path.root, recursive = T, showWarnings = F)

  xset.default <- do.call(xcmsSet, append(list(File_list), xcmsSet_param))

  xset.group <- xcms::group(xset.default, method="density", bw=bw_param[1], mzwid=mzwid_param, minfrac=minfrac_param)
  xset.2 <- xcms::retcor(xset.group, method="obiwarp", profStep=profStep_param, plottype="deviation")
  dev.copy(png, paste0(Results.path.root, "RetCor_01.png"), h=800, w=1600)
  graphics.off()
  xset.group2 <- xcms::group(xset.2, method="density", bw=bw_param[2], mzwid=mzwid_param, minfrac=minfrac_param)
  xset.3 <- xcms::retcor(xset.group2, method="obiwarp", profStep=profStep_param, plottype="deviation")
  dev.copy(png, paste0(Results.path.root, "RetCor_02.png"), h=800, w=1600)
  graphics.off()
  xset.group.4 <- xcms::group(xset.3, method="density", bw=bw_param[3], mzwid=mzwid_param, minfrac=minfrac_param)
  xset.filled <- xcms::fillPeaks(xset.group.4)

  ## Save peaks table
  Peak_Table_func <- peakTable(xset.filled)
  write.table(Peak_Table_func, file = paste0(Results.path.root, "Peak_Table.csv"), sep=";", col.names = NA)

  ## Increment table with parameters and results
  Parameters.Summary.temp <- data.frame(Groups = paste0(unique(xset.filled@phenoData$class), sep="", collapse = ", "),
                                        Sple.Nb = length(xset.filled@filepaths),
                                        Peak.Nb = nrow(xset.filled@peaks),
                                        Peak.Spl = round(nrow(xset.filled@peaks)/length(xset.filled@filepaths), 0),
                                        Pks.Grp.Nb = length(xcms.object@groupidx),
                                        Prof.Meth = xset.filled@profinfo[[1]],
                                        Prof.Step = xset.filled@profinfo[[2]],
                                        as.data.frame(t(unlist(xcmsSet_param[1:8]))))

  if (file.exists(paste0(Results.path, "Parameters.Summary.csv"))) {
    Parameters.Summary <- read.csv(file = paste0(Results.path, "Parameters.Summary.csv"), sep=";", row.names = 1)
    Parameters.Summary <- rbind(Parameters.Summary, Parameters.Summary.temp)
    rm(Parameters.Summary.temp)
    write.table(Parameters.Summary, file = paste0(Results.path, "Parameters.Summary.csv"), sep=";", col.names = NA)
  } else {
    Parameters.Summary <- Parameters.Summary.temp
    write.table(Parameters.Summary, file = paste0(Results.path, "Parameters.Summary.csv"), sep=";", col.names = NA)
    }

  if (QCs_Graph == TRUE) {
    png(filename = paste0(Results.path.root, "QCs.png"), h=1680, w=2400)
    par(mfrow=c(2,3))
    if(!is.null(xset.filled@phenoData$injectionOrder) & !is.null(xset.filled@phenoData$batch)) {
    Ordered_data <- xset.filled@phenoData[order(xset.filled@phenoData$injectionOrder),]
    plotQC(xset.filled, what="mzdevhist", sampNames = Ordered_data$injectionOrder, sampColors = Ordered_data$batch, sampOrder = order(xset.filled@phenoData$injectionOrder))
    plotQC(xset.filled, what="rtdevhist",   sampNames = Ordered_data$injectionOrder, sampColors = Ordered_data$batch, sampOrder = order(xset.filled@phenoData$injectionOrder))
    plotQC(xset.filled, what="mzdevmass",   sampNames = Ordered_data$injectionOrder, sampColors = Ordered_data$batch, sampOrder = order(xset.filled@phenoData$injectionOrder))
    plotQC(xset.filled, what="mzdevtime",   sampNames = Ordered_data$injectionOrder, sampColors = Ordered_data$batch, sampOrder = order(xset.filled@phenoData$injectionOrder))
    plotQC(xset.filled, what="mzdevsample", sampNames = Ordered_data$injectionOrder, sampColors = Ordered_data$batch, sampOrder = order(xset.filled@phenoData$injectionOrder))
    plotQC(xset.filled, what="rtdevsample", sampNames = Ordered_data$injectionOrder, sampColors = Ordered_data$batch, sampOrder = order(xset.filled@phenoData$injectionOrder))
    graphics.off()
    } else {
      Ordered_data <- xcms.object@phenoData[order(xcms.object@phenoData$injectionOrder),]
      plotQC(xset.filled, what="mzdevhist")
      plotQC(xset.filled, what="rtdevhist")
      plotQC(xset.filled, what="mzdevmass")
      plotQC(xset.filled, what="mzdevtime")
      plotQC(xset.filled, what="mzdevsample")
      plotQC(xset.filled, what="rtdevsample")
      graphics.off()
    }
  }

  if(is.data.frame(STDs_data) & !is.null(STDs_data$mz) & !is.null(STDs_data$ppm)) {

    Mz_ranges <- apply(STDs_data, 1, function(x) range(xcms:::ppmDev(as.numeric(x["mz"]), as.numeric(x["ppm"]))))
    for (i in ncol(Mz_ranges)){
      mz_range <- seq(Mz_ranges[1,i], Mz_ranges[2,i], by = 0.0001)
      temp <- subset(Peak_Table_func, round(mz, 4) %in% mz_range)
      if(exists("STD.subset")) { STD.subset <- rbind(STD.subset, temp) } else { STD.subset <- temp }
    }

    temp.plot <- melt(STD.subset, id.vars = c(1:(7+length(unique(xset.filled@phenoData$class)))))
    temp.plot2 <- merge(temp.plot, xset.filled@phenoData, by = "variable", all.x = T)
    write.table(STD.subset, file = paste0(Results.path.root, "Ions_Subset.csv"), sep=";", col.names = NA)
    write.table(temp.plot2, file = paste0(Results.path.root, "Ions_Subset_Metadata.csv"), sep=";", col.names = NA)

    temp_plot <- ggplot(temp.plot2, aes(x = as.factor(round(mz,2)), y = value, fill = variable, color = batch)) +
      geom_bar(stat="identity", position = "dodge") +
      ylab("") +
      xlab("") +
      ggtitle("") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      scale_fill_grey(start = 0, end = 1)

    png(filename = paste0(Results.path.root, "Ions_selection.png"), h=800, w=1600)
    print(temp_plot)
    graphics.off()

    if(STDs_EIC == TRUE) {
      nrow_val <- nrow(STD.subset)
      h <- 300*nrow_val
      png(filename = paste0(Results.path.root, "EIC_Ions_selection.png"), h=h, w=600)
      par(mfrow=c(nrow_val,1))
      plot(getEIC(xcms.object, groupidx = as.numeric(rownames(STD.subset)), rt = "corrected"), xcms.object)
      graphics.off()
    }

  } else { print("Need a dataframe with 'mz' column to analyse specifc ions") }



  return(xset.filled) ## Output results
}


#' xcms_orbi_A2
#'
#' This function perform the first steps of metabolomic analysis : xcmsSet and
#' two iteration with group and retcor to end with fillpeaks. The resulting output
#' is a xcms set object. Parameters are saved in a table.
#' @param File_list        File list to pass to xcmsSet method = 'centWave', ppm=7, peakwidth=c(4,20), snthresh=10, prefilter=c(4,1000), mzdiff=0.015, fitgauss=F, nSlaves = 4, phenoData = Sample.Metadata)
#' If phenoData is a dataframe with at lest a 'class' column. Could be any number of columns but in files sample order. If 'injectionorder' and 'batch" columns are
#' presents, those will be use to order en color QCs plot.
#' @param xcmsSet_param    xcmsSet parameters (can be any parameters formated like this : c(method="centWave", "ppm=7)).
#' @param Results.dir.name Name of the subfolder to store results.
#' @param bw_param         Vector for bw settings to use in 1, 2 and 3 iteration : c(1, 2, 3).
#' @param mzwid_param      mzwid parameter to use.
#' @param minfrac_param    minfrac parameter to use.
#' @param profStep_param   profStep parameter to use.
#' @param STDs_data        A dataframe of m/z to plot with at least one "mz" named column.
#' @param STDs_EIC         Should EIC of ions list be plotted.
#' @param QCs_Graph        Logical to determine of QCs graphs needs to be saved.
#' @keywords xcms, orbitrap
#' @usage xcms_orbi_A2()
#' xcms_orbi_A2()
#' @export

xcms_orbi_A2 <- function(File_list,
                        xcmsSet_param    = list(method="centWave",
                                                ppm=7,
                                                peakwidth=c(4,20),
                                                snthresh=10,
                                                prefilter=c(4,10000),
                                                mzdiff=-0.015,
                                                fitgauss=FALSE,
                                                nSlaves=4,
                                                phenoData = Sample.Metadata),
                        group_param      = list(method = "density",
                                                mzwid = 0.015,
                                                minfrac = 0.7),
                        group_bw         = c(15, 8, 0.8),
                        retcor_param     = list(method = "obiwarp",
                                                profStep = 1,
                                                plottype = "none"),
                        Results.dir.name = Results.path,
                        STDs_data        = NULL,
                        STDs_EIC         = FALSE,
                        QCs_Graph        = FALSE) {

  ## Package requirement
  require("xcms") ; require("reshape2") ; require("ggplot2")

  ## Create directory and path
  i <- 1
  Results.path.root <- paste0(Results.path, "XCMS_Result_", i, "/")
  while (dir.exists(Results.path.root)==TRUE) {
    i <- i+1
    Results.path.root <- paste0(Results.path, "XCMS_Result_", i, "/")
  }
  dir.create(path = Results.path.root, recursive = T, showWarnings = F)

  ## Workflow analysis
  xset.1 <- do.call(xcms::xcmsSet, append(list(File_list), xcmsSet_param))
  xset.group.1   <- do.call(xcms::group, append(append(alist(xset.1), group_param), list(bw = group_bw[1])))
  png(filename = paste0(Results.path.root, "RetCor.png"), h=1080, w=1080)
  par(mfrow=c(2,1))
  xset.2       <- do.call(xcms::retcor, append(alist(xset.group.1), retcor_param))
  xset.group.2 <- do.call(xcms::group, append(append(alist(xset.2), group_param), list(bw = group_bw[2])))
  xset.3       <- do.call(xcms::retcor, append(alist(xset.group.2), retcor_param))
  graphics.off()
  xset.group.3 <- do.call(xcms::group, append(append(alist(xset.3), group_param), list(bw = group_bw[3])))
  xset.filled  <- xcms::fillPeaks(xset.group.3)
  rm(xset.1, xset.group.1, xset.2, xset.group.2, xset.3, xset.group.3)

    ## Save peaks table
  Peak_Table_func <- peakTable(xset.filled)
  write.table(Peak_Table_func, file = paste0(Results.path.root, "Peak_Table.csv"), sep=";", col.names = NA)

  ## Increment table with parameters and results
  Parameters.Summary.temp <- data.frame("Groups" = paste0(unique(xset.filled@phenoData$class), sep="", collapse = ", "),
                                        "Sple.Nb" = length(xset.filled@filepaths),
                                        "Peak.Nb" = nrow(xset.filled@peaks),
                                        "Peak.Spl" = round(nrow(xset.filled@peaks)/length(xset.filled@filepaths), 0),
                                        "Pks.Grp.Nb" = length(xset.filled@groupidx),
                                        "Prof.Meth" = xset.filled@profinfo[[1]],
                                        "Prof.Step" = xset.filled@profinfo[[2]],
                                        "XcmsSet" = paste(unlist(xcmsSet_param), collapse = " "),
                                        "Group" = paste(unlist(group_param), collapse = " "),
                                        "Group_bw" = paste(unlist(group_bw), collapse = " "),
                                        "RetCor" = paste(unlist(retcor_param), collapse = " ")
                                        )

  ## Save parameters in dataframe
  if (file.exists(paste0(Results.path, "Parameters.Summary.csv"))) {
    Parameters.Summary <- read.csv(file = paste0(Results.path, "Parameters.Summary.csv"), sep=";", row.names = 1)
    Parameters.Summary <- rbind(Parameters.Summary, Parameters.Summary.temp)
    rm(Parameters.Summary.temp)
    write.table(Parameters.Summary, file = paste0(Results.path, "Parameters.Summary.csv"), sep=";", col.names = NA)
  } else {
    Parameters.Summary <- Parameters.Summary.temp
    write.table(Parameters.Summary, file = paste0(Results.path, "Parameters.Summary.csv"), sep=";", col.names = NA)
  }

  ## Generate QCs graph
  if (QCs_Graph == TRUE) {
    png(filename = paste0(Results.path.root, "QCs.png"), h=1680, w=2400)
    par(mfrow=c(2,3))
    if(!is.null(xset.filled@phenoData$injectionOrder) & !is.null(xset.filled@phenoData$batch)) {
      Ordered_data <- xset.filled@phenoData[order(xset.filled@phenoData$injectionOrder),]
      plotQC(xset.filled, what="mzdevhist", sampNames = Ordered_data$injectionOrder, sampColors = Ordered_data$batch, sampOrder = order(xset.filled@phenoData$injectionOrder))
      plotQC(xset.filled, what="rtdevhist",   sampNames = Ordered_data$injectionOrder, sampColors = Ordered_data$batch, sampOrder = order(xset.filled@phenoData$injectionOrder))
      plotQC(xset.filled, what="mzdevmass",   sampNames = Ordered_data$injectionOrder, sampColors = Ordered_data$batch, sampOrder = order(xset.filled@phenoData$injectionOrder))
      plotQC(xset.filled, what="mzdevtime",   sampNames = Ordered_data$injectionOrder, sampColors = Ordered_data$batch, sampOrder = order(xset.filled@phenoData$injectionOrder))
      plotQC(xset.filled, what="mzdevsample", sampNames = Ordered_data$injectionOrder, sampColors = Ordered_data$batch, sampOrder = order(xset.filled@phenoData$injectionOrder))
      plotQC(xset.filled, what="rtdevsample", sampNames = Ordered_data$injectionOrder, sampColors = Ordered_data$batch, sampOrder = order(xset.filled@phenoData$injectionOrder))
      graphics.off()
    } else {
      Ordered_data <- xset.filled@phenoData[order(xset.filled@phenoData$injectionOrder),]
      plotQC(xset.filled, what="mzdevhist")
      plotQC(xset.filled, what="rtdevhist")
      plotQC(xset.filled, what="mzdevmass")
      plotQC(xset.filled, what="mzdevtime")
      plotQC(xset.filled, what="mzdevsample")
      plotQC(xset.filled, what="rtdevsample")
      graphics.off()
    }
  }

  ## Filter data with specific mz
  if(is.data.frame(STDs_data) & is.numeric(STDs_data$mz) & is.numeric(STDs_data$ppm)) {
    Mz_ranges <- apply(STDs_data, 1, function(x) range(xcms:::ppmDev(as.numeric(x["mz"]), as.numeric(x["ppm"]))))
    for (i in 1:ncol(Mz_ranges)){
      temp <- subset(Peak_Table_func, mz > Mz_ranges[1,i] & mz < Mz_ranges[2,i])
      if(exists("STD.subset")) { STD.subset <- rbind(STD.subset, temp) } else { STD.subset <- temp }
    }

    if(nrow(STD.subset) > 0) {
      write.table(STD.subset, file = paste0(Results.path.root, "Ions_Subset.csv"), sep=";", col.names = NA)
      temp.plot <- melt(STD.subset, id.vars = c(1:(7+length(unique(xset.filled@phenoData$class)))))
      temp.plot2 <- merge(temp.plot, xset.filled@phenoData, by = "variable", all.x = T)
      write.table(temp.plot2, file = paste0(Results.path.root, "Ions_Subset_Metadata.csv"), sep=";", col.names = NA)

      temp_plot <- ggplot(temp.plot2, aes(x = as.factor(round(mz,2)), y = value, fill = variable, color = batch)) +
        geom_bar(stat="identity", position = "dodge") +
        ylab("") +
        xlab("") +
        ggtitle("") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
        scale_fill_grey(start = 0, end = 1)

      png(filename = paste0(Results.path.root, "Ions_selection.png"), h=800, w=1600)
      print(temp_plot)
      graphics.off()

      if(STDs_EIC == TRUE) {
        nrow_val <- nrow(STD.subset)
        h_param <- as.numeric(300 * nrow_val)
        png(filename = paste0(Results.path.root, "EIC_Ions_selection.png"), h = h_param, w=600)
        par(mfrow=c(nrow_val, 1))
        plot(getEIC(xset.filled, groupidx = as.numeric(rownames(STD.subset)), rt = "corrected"), xset.filled)
        graphics.off()
      }
    } else { print("No Specific ions found with given parameters.") }
  } else { print("Need a dataframe with 'mz' and 'ppm' column to analyse specifc ions") }

  return(xset.filled)
}







#' xcms_orbi_Results
#'
#' This function generate peak_table, STDs EIC across samples, PCA (optional) and return a list with [1] datamatrix, [2] sample.metadata and [3] variable.metadata.
#' @param filled_peak_object An xcmsSet object with filled peaks
#' @param STDs_mass Vector of STDs exact mass
#' @param STDs_ppm  Deviation for STDs in ppm (try to increase if your STDs are not found)
#' @param perform_PCA Logical to do PCA analysis or not
#' @param Sample.Metadata Dataframe with samples metadata
#' @param PCA_group vector with samples metadata column to use for groups in ACP
#' @keywords xcms, orbitrap
#' @usage xcms_orbi_GRT(xcms_set_obj, Results.dir.name="Default", bw_param=c(25, 10, 0.7), mzwid_param=0.005, minfrac_param=0.25, profStep_param=0.8)
#' xcms_orbi_GRT()
#' @export

xcms_orbi_Results <- function(filled_peak_object,
                              Results.dir.name = "Default",
                              STD = c(TRUE, FALSE),
                              STDs_mass = c(133.1062, 206.1014, 179.0874, 281.3265),
                              STDs_ppm = 10,
                              perform_PCA = c(TRUE, FALSE),
                              Sample.Metadata,
                              PCA_group = c(1,2,3,4)){
  ## Package requirement
  library("xcms")

  ## Create directory and path
  Results.path.root <- paste0("./", Results.dir.name, "/")
  dir.create(Results.path.root, showWarnings = F, recursive = T)

  ## Get group metadata and sample data
  write.table(peakTable(filled_peak_object), file=paste0(Results.path.root, "Peak_Table.csv"), sep=";", col.names=NA)
  ## Generate files for opls analysis
  Data <- list()
  Peak_Table_func.temp <- peakTable(filled_peak_object)
  Data$Datamatrix <- data.frame(t(Peak_Table_func.temp[(ncol(Peak_Table_func.temp)-nrow(filled_peak_object@phenoData)+1):ncol(Peak_Table_func.temp)]))
  if (is.data.frame(Sample.Metadata)==T){
    Data$Sample.metadata <- subset(Sample.Metadata, rownames(Sample.Metadata) %in% rownames(filled_peak_object@phenoData))
  } else {
    print("Sample.Metadata is not a data.frame, return empty list[[2]]")
    Data$Sample.metadata <- NULL
  }
  Data$Variable.metadata <- Peak_Table_func.temp[1:(ncol(Peak_Table_func.temp)-nrow(filled_peak_object@phenoData))]

  ## Extract Standards infos
  if (STD == TRUE) {
    Results.path.std <- paste0(Results.path.root, "STD/")
    dir.create(Results.path.std, showWarnings = F)
    require(xcms)
    group.id.p <- c()
    for (i in 1:length(STDs_mass)){
      temp <- rownames(subset(Data[[3]], Data[[3]][, "mz"] >= min(xcms:::ppmDev(STDs_mass[i], STDs_ppm)) & Data[[3]][, "mz"] <= max(xcms:::ppmDev(STDs_mass[i], STDs_ppm))))
      group.id.p <- union(group.id.p, temp)
    }
    group.id.p <- as.numeric(group.id.p)
    temp.eic.r <- getEIC(filled_peak_object, groupidx=group.id.p, rt="raw")
    temp.eic.c <- getEIC(filled_peak_object, groupidx=group.id.p, rt="corrected")
    png(filename=paste0(Results.path.std, "EIC.STD.peaks_", ".png"), width=800, height=1080, units="px", res=100)
    par(mfrow=c(length(group.id.p),2))
    for (i in 1:length(group.id.p)){
      plot(temp.eic.r, filled_peak_object, groupidx=i, main="RAW")
      plot(temp.eic.c, filled_peak_object, groupidx=i, main="Corrected")
    }
    dev.off()

    STDs.results <- subset(Data[[3]], rownames(Data[[3]]) %in% group.id.p)
    STDs.results$ppm <- (STDs.results$mzmax - STDs.results$mzmin)/(STDs.results$mz/1000000)
    write.table(STDs.results, file=paste0(Results.path.std, "Results_STDs.csv"), sep=";", col.names = NA)
    print(STDs.results)
  }

  ## Perform PCA
  if (perform_PCA == TRUE) {
    Results.path.pca <- paste0(Results.path.root, "PCA/")
    dir.create(Results.path.pca, showWarnings = F)
    library("ropls")
    ACP.results  <- opls(Data[[1]], predI=NA, plotL=F)
    png(filename=paste0(Results.path.pca,"ACP summary.png"), width=800, height=1200, units="px", res=150)
    par(mfrow=c(3,2))
    plot(ACP.results, typeVc="overview", parDevNewL=F)
    plot(ACP.results, typeVc="x-loading", parDevNewL=F)
    plot(ACP.results, typeVc="x-score", parDevNewL=F)
    plot(ACP.results, typeVc="outlier", parDevNewL=F)
    plot(ACP.results, typeVc="correlation", parDevNewL=F)
    dev.off()

    if (!is.null(PCA_group)) {
      if(length(names(Sample.Metadata)) >= length(PCA_group)) {
      sqrt.group <- sqrt(length(PCA_group))
      if (sqrt.group==round(sqrt.group)){
        x <- sqrt.group
        y <- sqrt.group
      } else {
        x <- round(sqrt.group,0)
        y <- ceiling(sqrt.group)
      }

      png(filename=paste0(Results.path.pca, "ACP_Ellipses.png"), width=x*300, height=y*300, units="px", res=100)
      par(mfrow=c(x, y))
      for (i in PCA_group){
        temp.factor <- Data[[2]][, i]
        temp.factor.names <- names(Data[[2]][i])
        plot(ACP.results, typeVc="x-score", parAsColFcVn=addNA(as.factor(temp.factor)), parEllipses=F, parDevNewL=F)
        text(par()$usr[1]/1.2, par()$usr[3]/1.1, temp.factor.names)
      }
      dev.off()
      }
    }
  }
  return(Data)
}


#' SD_batch_list
#'
#' This function separate .mzXML file list by batch. Useful when you use directory for your classes
#' and want to do separate xcms analysis by another grouping factor like batch sequences.
#' @param Files.dir Path to your results from workdirectory.
#' @param Batch.list A dataframe with rownames is files names (without extension) and first column is batch number (or other factor).
#' @keywords xcms, orbitrap
#' @export
#'
SD_batch_list <- function(Files.dir="./",
                          Batch.list=Batch.list){
  Files <- list.files(Files.dir, recursive=T, full.names=T)
  result.list <- list()
  for (i in unique(Batch.list[,1])){
    result.list[[i]] <- subset(Files, gsub("*.mzXML","", x=basename(Files)) %in% rownames(subset(Batch.list, Batch.list[1]==i)))
  }
  return(result.list)
}


#' SD_batch_set
#'
#' This function separate .mzXML file list by batch. Useful when you use directory for your classes
#' and want to do separate xcms analysis by another grouping factor like batch sequences.
#' @param Batch.files List returned by SD_mass_batch with files for each batch.
#' @param xcmsSet_param xcmsSet parameters (can be any parameters formated like this : c(method="centWave", "ppm=7))
#' @keywords xcms, orbitrap
#' @export
#'
SD_batch_set <- function(Batch.files=Batch.files,
                         xcmsSet_param = list(method="centWave", ppm=7, peakwidth=c(4,20), snthresh=10, prefilter=c(4,10000), mzdiff=-0.001, fitgauss=FALSE, nSlaves=4)
                         ){
  Batch.xcmset.list <- list()
  for (i in 1:length(Batch.files)){
    a <- paste0("Batch.", i)
    Batch.xcmset.list[[i]] <- SDjoygret::xcms_orbi_A(Batch.files[[i]], Results.dir.name = a, xcmsSet_param = xcmsSet_param)
  }
  return(Batch.xcmset.list)
}


#' SD_files_class
#'
#' This function get multiple files path and name, and then parse the information into a dataframe,
#' making files grouping by directory path easy.
#' @param dir_path Origin directory pathway.
#' @keywords Files, info
#' @examples
#' SD_files_class(dir_path = "./")
#' @export
#'
SD_files_class <- function(dir_path) {
  Files <- list.files("dir_path", recursive=T, full.names=T)
  Files.str <- as.data.frame(read.table(textConnection(sub("*.mzXML","", x=Files)), sep="/"))
  return(Files.str)
}


#' SD_pca
#'
#' Perform pca according to ropls method and save results in a folder.
#' You can choose to color samples by factor or value contained in a sample metadata table.
#' @param Data list with [[1]] Datamatrix [[2]] samples metadata
#' @param Results.path.root Name for the results folder
#' @param ropls_param Parameters to pass to ropls::opls function
#' @param factor_group vector of column(s) number/name of metadata table (in Data[[2]]) to use for grouping corcircles
#' @keywords pca
#' @export
#' @examples
#' SD_pca(Data, Results.path.root = "Default", ropls_param = list(predI = 2, plotL = F), factor_group = NULL)

SD_pca <- function(Data_pca,
                   Folder_name = "Default",
                   ropls_param = list(predI = 2, plotL = F),
                   factor_group = NULL
                   )
{
  Results.path <- paste0("./", Folder_name, "/")
  dir.create(Results.path, showWarnings = F)
  require(ropls)
  pca_result <- do.call(ropls::opls, append(list(x=(Data[[1]])), ropls_param))

  png(filename = paste0(Results.path, "ACP_result.png"), height = 3 * 400, width = 2 * 400, units = "px", res = 100)
  par(mfrow=c(3,2))
  plot(pca_result, typeVc = "overview", parDevNewL = F)
  plot(pca_result, typeVc = "x-loading", parDevNewL = F)
  plot(pca_result, typeVc = "x-score", parDevNewL = F)
  plot(pca_result, typeVc = "outlier", parDevNewL = F)
  plot(pca_result, typeVc = "correlation", parDevNewL = F)
  dev.off()
  if (!is.null(factor_group) & is.data.frame(Data_pca[[2]])){
    for (i in factor_group){
      temp.factor <- Data_pca[[2]][,i]
      temp.factor.names <- names(Data_pca[[2]][i])
      png(filename=paste0(Results.path, "ACP_Ellipses_", temp.factor.names, ".png"), height = 800, width = 800, units = "px", res = 100)
      plot(pca_result, typeVc = "x-score", parAsColFcVn=addNA(as.factor(temp.factor)), parEllipses = F, parDevNewL = F)
      text(par()$usr[1]/1.2, par()$usr[3]/1.1, temp.factor.names)
      dev.off()
    }
  }
  return(invisible(pca_result))
}


#' SD_subset_zero
#'
#' Subset a dataframe containing numeric value by selecting column which have less than
#' Z.Seuil percents of zero values.
#' @param x Dataframe (a datamatrix) with only numeric values
#' @param Z.Seuil Percentage threshold of zero values to accept for each variable
#' @keywords pca
#' @export
#' @examples
#' SD_subset_zero(Data[[1]], Z.Seuil = 50)

SD_subset_zero <- function(x, Z.Seuil = 80) {
  require(ggplot2)

  ## Calculate percentage of zero value in each column
  Zero.perc <- apply(x, 2, function(x) round(length(which(x == 0)) * 100 / length(x), 3))
  Zero.perc.2 <- data.frame(Zero.perc)

  ## Subset the selected column in the datamatrix
  Zer.perc.var <- subset(data.frame(Zero.perc), Zero.perc <= Z.Seuil)
  data.subset <- subset(x, select = rownames(Zer.perc.var))

  ## Plot graph of results and selected thresold
  Variable.Zero.Percentage.Plot <- ggplot(Zero.perc.2, aes(x = reorder(rownames(Zero.perc.2), Zero.perc))) +
    geom_point(aes(y = Zero.perc, color = rownames(Zero.perc.2) %in% rownames(Zer.perc.var)), size = 3, shape=21, fill = "white") +
    geom_hline(yintercept = Z.Seuil, color = "black", alpha = 0.5) +
    ylim(0,100) +
    theme(legend.position = 0) +
    coord_flip()

  ## Print results
  print(paste0("Variable(s) selected : ", length(rownames(Zer.perc.var)), " on ", length(rownames(x))))
  print(summary(Data.raw.zero.subset))

  return(invisible(list("Data.Subset" = data.subset,
                        "Zer.Perc" = Zero.perc,
                        "Plot" = Variable.Zero.Percentage.Plot,
                        "Zero.Threshold" = Z.Seuil)))
}





#' SD_data_sub
#'
#' Subset a list of three dataframes : [[1]] Datamatrix, [[2]] Samples.Metadata [[3]] Variable.Metadata.
#' [[1]] With samples ids as rownames and variables ids as column names
#' [[2]] Metadata of sample in column and samples ids as rownames (same order than rows in [[1]])
#' [[3]] Variables metadata in column and variable ids as rownames (same order than columns in [[1]])
#' @param Data list with [[1]] Datamatrix [[2]] samples metadata
#' @param Results.path.root Name for the results folder
#' @param ropls_param Parameters to pass to ropls::opls function
#' @param factor_group vector of column(s) number/name of metadata table (in Data[[2]]) to use for grouping corcircles
#' @keywords pca
#' @examples
#' SD_pca(Data, Results.path.root = "Default", ropls_param = list(predI = 2, plotL = F), factor_group = NULL)

Data_sub <- function(Data, col_select, ...)
{
  ## Check entry file
  if (missing(Data)) stop("arg. Data is missing.")
  if (!is.list(Data)) stop("arg. Data must be a list of two to three dataframe")
  if (length(Data) < 2) stop("arg. Data must be a list of two to three dataframe")
  if (length(Data) > 3) stop("arg. Data must be a list of two to three dataframe")
  invisible(lapply(Data, function(x) if(!is.data.frame(x)) stop("arg. Data must be a list of dataframe")))
  invisible(lapply(Data[[1]], function(x) if(!is.numeric(x)) stop("Data[[1]] need to be a matrix with numeric value")))

  ## Filter by sample metadata
  temp.choice <- c("A", "B")
  col_name <-

  Data.subset <- list()
  Data.subset[[2]] <- subset(Data[[2]], colnames(Data[[2]]) %in% temp.choice)
  Data.subset[[1]] <- subset(Data[[1]], rownames(Data[[1]]) %in% rownames(Data.subset[[2]]))
  ## Filter by variable metadata
  ## Filter by datamatrix values
}


#' get_sign
#'
#' Get the theoritical significance glm model.
#' This function come from WinVector (Nina Zumel, http://winvector.github.io/VariableSignal/VariableSignal.html).
#' @param a model generated with glm (use of $null.deviance, $deviance, df.null and df.residual)
#' @keywords glm, significance
#' @examples
#' get_sign(model)

get_sign = function(model) {
  delta_deviance = model$null.deviance - model$deviance
  df = model$df.null - model$df.residual
  sig = pchisq(delta_deviance, df, lower.tail=FALSE)
}






#' Import Excel file to list
#'
#' Import excel spreadsheet to a 3 levels list. This function doesn't work for big files because of memory issue with Java.
#' @param filename Path to the excel file to import
#' @return A 3 levels list
#' @keywords list, import, excel
#' @export
#' @examples
#' importWorksheets()

importWorksheets <- function(filename) {
  ## Sources : http://stackoverflow.com/questions/12945687/how-to-read-all-worksheets-in-an-excel-workbook-into-an-r-list-with-data-frame-e
  require(XLConnect)
  # filename: name of Excel file
  workbook <- loadWorkbook(filename)
  setMissingValue(workbook,  value = c("NA"))
  sheet_names <- getSheets(workbook)
  names(sheet_names) <- sheet_names
  sheet_list <- lapply(sheet_names, function(.sheet){
    readWorksheet(object=workbook, .sheet)})
  sheet_list <- lapply(sheet_list, function(x) {
    rownames(x) <- x[,1]
    x <- x[-1]
  })
  return(sheet_list)
}




#' Check format of 3 levels list
#'
#' Check the format of 3 levels list.
#' @param data Three levels list with [[1]] Datamatrix, [[2]] SamplesMetadata, [[3]] VariableMetadata.
#' @keywords list, check
#' @export
#' @examples
#' check.list.format()

check.list.format <- function (data) {
  if(!is.list(List.Result)){stop("Data should be a list with (1) Datamatrix (2) Sample.Metadata (3) Variable.Metadata")}
  if(FALSE %in% c(lapply(List.Result, class) == "data.frame")){stop("List levels should be data.frame")}
  if(!identical(colnames(List.Result[[1]]), rownames(List.Result[[3]]))){stop("Datamatrix colnames should be identical of Variable.Metadata rownames")}
  if(!identical(rownames(List.Result[[1]]), rownames(List.Result[[2]]))){stop("Datamatrix rownames should be identical of Sample.Metadata rownames")}
  dim.temp <- lapply(List.Result, dim)
  if(!dim.temp[[1]][1] == dim.temp[[2]][1]){stop("Datamatrix row number should be the same as Sample.Metadata")}
  if(!dim.temp[[1]][2] == dim.temp[[3]][1]){stop("Datamatrix col number should be the same as Variable.Metadata row number")}
  print("Data seems OK")
}


#' Subset list
#'
#' Subset a list of three dataframes : [[1]] Datamatrix, [[2]] SamplesMetadata [[3]] VariableMetadata.
#' [[1]] samples.ID x variables.ID
#' [[2]] samples.ID x sample.metadata (same order than rows in [[1]])
#' [[3]] variables.ID x variable.metadata (same order than columns in [[1]])
#' @param data list of three dataframes : [[1]] Datamatrix, [[2]] SamplesMetadata [[3]] VariableMetadata.
#' @param Var.sel vector of variable to subset (rownames) from [[3]]
#' @param Samples.sel vector of samples to subset (rownames) from [[2]]
#' @keywords subset, list
#' @export
#' @examples
#' data.subset()
#'

data.subset <- function(data, Var.sel = NULL, Samples.sel = NULL) {
  check.list.format(data)
  temp.data <- data
  if(!is.null(Samples.sel)){
    temp.data[[2]] <- subset(temp.data[[2]], rownames(temp.data[[2]]) %in% Samples.sel)}
  if(!is.null(Var.sel)){
    temp.data[[3]] <- subset(temp.data[[3]], rownames(temp.data[[3]]) %in% Var.sel)}
  if(length(data)==3){
    temp.data[[1]] <- subset(temp.data[[1]], rownames(temp.data[[1]]) %in% rownames(temp.data[[2]]), select = rownames(temp.data[[3]]))
    return(temp.data)} else {
      warning("No Variable metadata")
      temp.data[[1]] <- subset(temp.data[[1]], rownames(temp.data[[1]]) %in% rownames(temp.data[[2]]))
      return(temp.data)}
}


#' Split a dataframe
#'
#' Split a dataframe in two.
#' @param dataframe a dataframe to subset.
#' @param p the proportion of data to keep in trainset
#' @param seed the seed number for repetability (same seed will generate same subset for a given dataframe).
#' @keywords dataframe, training, subset
#' @return A list with [1] trainset dataframe, [2] testset dataframe.
#' @usage splitdf(cars, p = 0.5, seed = 95687)
#' @export
#' @examples
#' splitdf()

splitdf <- function(dataframe, p = 0.5, seed = 95687) {
  if (!is.null(seed)) set.seed(seed)
  index <- 1:nrow(dataframe)
  trainindex <- sample(index, trunc(length(index)*p))
  trainset <- dataframe[trainindex, ]
  testset <- dataframe[-trainindex, ]
  return(list(trainset=trainset,testset=testset))
}


#' Write CSV for Excel Fr
#'
#' Write a csv directly sees as a table by excel.
#' @param data Dataframe to save
#' @param file File path and name
#' @keywords dataframe, excel, save
#' @return Write a csv file
#' @usage write.csv3(data, file = "./Results/Table_x.csv")
#' @export
#' @examples
#' write.csv3()
write.csv3 <- function(data, file) {
  write.table(data, file = file,
              sep = ";",
              dec = ".",
              append = F,
              qmethod = "double",
              row.names = FALSE)
}


#' Save file for GALAXY
#'
#' Write 3 csv files from a three levels list to import in GALAXY.
#' @param x 3 levels list
#' @param Results.path Saving path
#' @param pref Prefix to add
#' @keywords list, galaxy, export
#' @return Write 3 csv file on disk
#' @usage galaxy.save.list(data.list)
#' @export
#' @examples
#' galaxy.save.list()
galaxy.save.list <- function(x, Results.path = "./", pref = "GALAXY-"){
  List.format.check(x)
  temp <- t(x[[1]])
  write.table(data.frame("dataMatrix" = rownames(temp), temp), file = paste0(Results.path, pref, "Datamatrix.csv"), sep = "\t", quote = F, row.names = F)
  write.table(data.frame("sampleMetadata" = rownames(x[[2]]), x[[2]]), file = paste0(Results.path, pref, "SampleMetadata.csv"), sep = "\t", quote = F, row.names = F)
  write.table(data.frame("variableMetadata" = rownames(x[[3]]), x[[3]]), file = paste0(Results.path, pref, "VariableMetadata.csv"), sep = "\t", quote = F, row.names = F)
}


#' Get legend of a ggplot
#'
#' Get the legend of a ggplot to draw it externaly & draw it as grob.
#' @param a.gplot A ggplot
#' @keywords legend, ggplot
#' @return Write 3 csv file on disk
#' @usage galaxy.save.list(data.list)
#' @export
#' @examples
#' g_legend()
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}


#' Find scales for plots
#'
#' Return the limits of two vectors centered on 0 and with a 10 % margin.
#' @param x Values on the x axis
#' @param y Values on the y axis
#' @keywords scales, ggplot
#' @return 4 numerics with max and min value for x and y
#' @usage find.limits(x,y)
#' @export
#' @examples
#' find.limits()
find.limits <- function(x,y) {
  temp.list <- rbind(c(max(abs(x))*-1.1, max(abs(x))*1.1),
                     c(max(abs(y))*-1.1, max(abs(y))*1.1))
  return(temp.list)
}


#' Semi-automatic ggplot theme
#'
#' Returns the ploting parameters according to input, with working default for missing ones.
#' Used for custom ggplot of multivariate results (pca, pls, opls)
#' @param Samples.grp Name of the grouping factor for samples
#' @param Variables.grp Name of the grouping factor for variables
#' @param limits Limits of the axis (returned by find.limits)
#' @param Legend.L Logical to draw legend
#' @param colorL Logical to use color or grey scale
#' @param labels list for titles (title, x, y)
#' @param geom_path Logical for drawing path
#' @param labelsL Add labels to points
#' @keywords scales, ggplot
#' @return ggplot list of aestethic
#' @export
#' @examples
#' plot.theme.1()
plot.theme.1 <- function(Samples.grp = NULL, Variables.grp = NULL, limits = NULL, Legend.L = T, colorL = F, labels = list(title = "", x = "", y = ""), geom_path = F, labelsL = F) {
  opls.ggplot_theme <- list(
    if(!is.null(Samples.grp) & geom_path == T){geom_path(alpha = 0.4)},
    if(!is.null(Samples.grp)){labs(colour = Samples.grp)},
    if(!is.null(Samples.grp)){aes(group = as.factor(get(Samples.grp)), color = as.factor(get(Samples.grp)))},
    if(!is.null(Variables.grp)){aes(color = as.factor(get(Variables.grp)))},
    if(!is.null(Variables.grp)){labs(color = Variables.grp)},
    if(colorL == F) {scale_colour_grey()},
    geom_hline(yintercept = 0, linetype = 2, color = "grey"),
    geom_vline(xintercept = 0, linetype = 2, color = "grey"),
    geom_point(alpha = 0.8),
    if(!is.null(limits)){xlim(limits[1,1], limits[1,2])},
    if(!is.null(limits)){ylim(limits[2,1], limits[2,2])},
    ggplot_theme,
    labs(labels),
    if(labelsL == T){geom_text(vjust = -0.8)},
    if(Legend.L == F){theme(legend.position = 0)} else {theme(legend.position = c(0,0), legend.justification = c(0,0), legend.direction = "horizontal", legend.title = element_blank())}
  )
  return(opls.ggplot_theme)
}


#' Perform PCA and custom plot
#'
#' Perform a PCA analysis on a 3 levels list and create custom plot.
#' @param Data.list A 3 levels list with the data
#' @param Samples.grp Factor name used for grouping samples (affect plot only)
#' @param Variables.grp Factor name used for grouping variables (affect plot only)
#' @param Legend.L Logical for drawing legends
#' @param colorL Logical for using color or greyscale
#' @param Samp.lab.L Logical for drawing Sample labels
#' @param Var.lab.L Logical for drawing Variable labels
#' @keywords pca, ggplot
#' @return A list with [1] pca results, [2] plot as grobs.
#' @export
#' @examples
#' datamatrix.pca()
datamatrix.pca <- function (Data.list, Samples.grp = NULL, Variables.grp = NULL, Legend.L = F, colorL = F, Samp.lab.L = F, Var.lab.L = T) {
  require(ropls) ; require(ggplot2) ; require(gridExtra)
  temp.pca <- opls(Data.list[[1]], predI = 2, plotL = F)
  temp.scores <- merge(Data.list[[2]], data.frame(temp.pca$scoreMN), by.x = 0, by.y = 0)
  limits1 <- find.limits(temp.scores$p1,temp.scores$p2)
  temp.loadings <- merge(Data.list[[3]], data.frame(temp.pca$loadingMN), by.x = 0, by.y = 0)
  limits2 <- find.limits(temp.loadings$p1,temp.loadings$p2)
  ## plot variables
  labels1 <- list(title = paste0("Scores plot ", temp.pca$descriptionMC[1], " samples\n(", temp.pca$descriptionMC[4], " missing values)"),
                  x = paste0("p1 (", temp.pca$modelDF$R2X[1]*100, " %)"),
                  y = paste0("p2 (", temp.pca$modelDF$R2X[2]*100, " %)"))
  labels2 <- list(title = paste0("Loadings plot\n", temp.pca$descriptionMC[2], " variables (", temp.pca$descriptionMC[3], " excluded)"),
                  x = paste0("p1 (", temp.pca$modelDF$R2X[1]*100, " %)"),
                  y = paste0("p2 (", temp.pca$modelDF$R2X[2]*100, " %)"))
  ## plots
  plot1 <- ggplot(temp.scores, aes(p1, p2)) +
    plot.theme.1(Samples.grp, Variables.grp = NULL, limits = limits1, Legend.L, colorL, labels1, geom_path = T, labelsL = Samp.lab.L)
  plot2 <- ggplot(temp.loadings, aes(p1, p2, label = Row.names)) +
    plot.theme.1(Samples.grp = NULL, Variables.grp, limits = limits2, Legend.L, colorL, labels1, geom_path = F, labelsL = Var.lab.L)
  ## Result
  return(list("PCA" = temp.pca, "Plot" = grid.arrange(plot1, plot2, nrow = 1 )))
}









