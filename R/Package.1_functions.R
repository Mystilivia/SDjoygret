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
  require("grDevices")
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
  require("grDevices")
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
  require("xcms") ; require("reshape2") ; require("ggplot2") ; require("grDevices")
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
  library("xcms") ; require("grDevices")

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


#' Check format of dlist
#'
#' Check the format of 3 levels list.
#'
#' Chech if the list elements are data.table or data.frame, id rownames [[1]] and [[2]] are identical (for
#' data.frame only, use first column as row id with data.table), and if colnames [[1]] are identical to
#' rownames [[3]]. Check also dataframes dimension to check consistency and propose to transform data.frame
#' to data.table.
#'
#' @param dlist Three levels list with [[1]] Datamatrix, [[2]] SamplesMetadata, [[3]] VariableMetadata.
#' @param to.data.table Logical to convert data.frame's dlist to data.table with to.data.table function
#' @return Print result of check as character and return the dlist (or converted dlist) if format is ok
#' @keywords list, check
#' @export
#' @examples
#' check.list.format()
check.list.format <- function (dlist, to.data.table.L = T) {
  require("data.table")
  if(!is.list(dlist)){stop("Data should be a list with (1) Datamatrix (2) Sample.Metadata (3) Variable.Metadata")}
  temp.data.str <- data.table("ListLevel" = names(dlist),
                              "matrix" = lapply(dlist, function(x) {any(class(x) == "matrix")}),
                              "data.table" = lapply(dlist, function(x) {any(class(x) == "data.table")}),
                              "data.frame" = lapply(dlist, function(x) {any(class(x) == "data.frame")}),
                              "tibble" = lapply(dlist, function(x) {any(class(x) == "tibble")}))
  if(temp.data.str[,any(data.table == F)] & temp.data.str[,any(data.frame==F)]) {stop("List levels should be data.frame or data.table")}
  if(!temp.data.str[,any(data.table == F)]) { ## all are data.table
    if(!identical(names(dlist[[1]])[-1], dlist[[3]][[1]])){stop("Datamatrix colnames should be identical of Variable.Metadata rownames")}
    if(!identical(dlist[[1]][[1]], dlist[[2]][[1]])){stop("Datamatrix rownames should be identical of Sample.Metadata rownames")}
    dim.temp <- lapply(dlist, dim)
    if(!dim.temp[[1]][1] == dim.temp[[2]][1]){stop("Datamatrix row number should be the same as Sample.Metadata")}
    if(!dim.temp[[1]][2]-1 == dim.temp[[3]][1]){stop("Datamatrix col number should be the same as Variable.Metadata row number")}
    print("Data format is ok")
    print(temp.data.str)
    return(dlist)
  } else {
    if(!temp.data.str[,any(data.frame==F)]) { ## all are data.table
      if(!identical(colnames(dlist[[1]]), rownames(dlist[[3]]))){stop("Datamatrix colnames should be identical of Variable.Metadata rownames")}
      if(!identical(rownames(dlist[[1]]), rownames(dlist[[2]]))){stop("Datamatrix rownames should be identical of Sample.Metadata rownames")}
      dim.temp <- lapply(dlist, dim)
      if(!dim.temp[[1]][1] == dim.temp[[2]][1]){stop("Datamatrix row number should be the same as Sample.Metadata")}
      if(!dim.temp[[1]][2] == dim.temp[[3]][1]){stop("Datamatrix col number should be the same as Variable.Metadata row number")}
      print("Data seems OK but are stored in data.frame")
      if(to.data.table.L == T) {
        to.data.table(dlist, rownames = T)
      } else {
        print("Data weren't converted to data.table, set arg to.data.table = T if convertion is needed")
        print(temp.data.str)
        return(dlist)
      }
    }
  }
}


#' Check dlist table format and convert to data.table
#'
#' Check the format of 3 levels list tables print a class summary and convert to
#' data.table if not. User can specify if there is rownames (if data.frame were used in entry) so
#' the function can place them as first variable (since data.table does not use rownames logic).
#'
#' @param dlist Three levels list with [[1]] Datamatrix, [[2]] SamplesMetadata, [[3]] VariableMetadata.
#' @param rownames Logical to specify wether the input dlist use data.frame rownames to store rowID.
#' @return The dlist with data converted to dataframe. Also print a check status of class used in the
#' dlist
#' @keywords list, check, data.table
#' @export
#' @examples
#' to.data.table()
to.data.table <- function(dlist, rownames = F) {
  require(data.table)
  temp.data.str <- data.table("ListLevel" = names(dlist),
                              "matrix" = lapply(dlist, function(x) {any(class(x) == "matrix")}),
                              "data.table" = lapply(dlist, function(x) {any(class(x) == "data.table")}),
                              "data.frame" = lapply(dlist, function(x) {any(class(x) == "data.frame")}),
                              "tibble" = lapply(dlist, function(x) {any(class(x) == "tibble")}))
  if(!temp.data.str[,any(data.table==F)]) { ## all are data.table
    print("Data seems ok")
    print(temp.data.str)
    return(dlist)
  } else {
    if(!temp.data.str[,any(data.frame==F)]) { ## at least one isn't a data.table but all are data.frame
      if(rownames == F) {
        Data.2 <- lapply(dlist, function(x) {
          data.table(x)
        })
        print("Data were converted to data.table")
        to.data.table(Data.2)
      } else {
        Data.2 <- lapply(dlist, function(x) {
          data.table("RowID" = rownames(x), x)
        })
        print("Data were converted to data.table and rownames added in RowID")
        to.data.table(Data.2, rownames = F)
      }
    } else {
      print("Data were not data.frame or data.table :")
      print(temp.data.str)
      return(NULL)
    }
  }
}




#' Import Excel file to list
#'
#' Import each sheet of an excel file to a list. This function use the readxl package, which seems
#' to handle large files. By default, NA is set to "NA" first column to rownames and first
#' row to column names. Rownames can be override with fcolL argument.
#' @param Data.path Path to the excel file to import
#' @param fcolL Logical to specify if th first column should be use for rownames (must contains unique names)
#' @return A list of dataframe
#' @keywords list, import, excel
#' @export
#' @examples
#' importWorksheets.xls()
importWorksheets.xls <- function(Data.path, fcolL = F) {
  require ("readxl")
  Sheet.names <- excel_sheets(Data.path)
  data.list <- list()
  for (i in Sheet.names) {
    data.list[[i]] <- data.frame(read_excel(Data.path, sheet = i, col_names = T, na = "NA"))
    if(isTRUE(fcolL)) {
    rownames(data.list[[i]]) <- data.list[[i]][[1]]
    data.list[[i]] <- data.list[[i]][-1]
    }
  }
  return(data.list)
}


#' Subset list
#'
#' Subset a three level list by variables and/or samples.
#' @param dlist list of three dataframes : [[1]] Datamatrix, [[2]] SamplesMetadata [[3]] VariableMetadata.
#' @param Var.sel vector of variable to subset (rownames) from [[3]]
#' @param Samples.sel vector of samples to subset (rownames) from [[2]]
#' @keywords subset, list
#' @export
#' @examples
#' dlist.subset()
dlist.subset <- function(dlist,
                         Var.sel = NULL,
                         Samples.sel = NULL) {
  check.list.format(dlist)
  temp.dlist <- dlist
  if(!is.null(Samples.sel)){
    temp.dlist[[2]] <- subset(temp.dlist[[2]], rownames(temp.dlist[[2]]) %in% Samples.sel)}
  if(!is.null(Var.sel)){
    temp.dlist[[3]] <- subset(temp.dlist[[3]], rownames(temp.dlist[[3]]) %in% Var.sel)}
  if(length(dlist)==3){
    temp.dlist[[1]] <- subset(temp.dlist[[1]], rownames(temp.dlist[[1]]) %in% rownames(temp.dlist[[2]]), select = rownames(temp.dlist[[3]]))
    return(temp.dlist)} else {
      warning("No Variable metadlist")
      temp.dlist[[1]] <- subset(temp.dlist[[1]], rownames(temp.dlist[[1]]) %in% rownames(temp.dlist[[2]]))
      return(temp.dlist)}
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
splitdf <- function(dataframe,
                    p = 0.5,
                    seed = 95687) {
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
write.csv3 <- function(data,
                       file) {
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
#' @param Results.path Saving path
#' @param pref Prefix to add
#' @inheritParams dlist.subset
#' @keywords list, galaxy, export
#' @return Write 3 csv file on disk
#' @usage galaxy.save.list(data.list, Results.path = "./", pref = "GALAdlistY-")
#' @export
#' @examples
#' galaxylist.save.list()
galaxylist.save.list <- function(dlist,
                             Results.path = "./",
                             pref = "GALAdlistY-"){
  check.list.format(dlist)
  temp <- t(dlist[[1]])
  write.table(data.frame("Datamatrix" = rownames(temp), temp), file = paste0(Results.path, pref, "Datamatridlist.csv"), sep = "\t", quote = F, row.names = F)
  write.table(data.frame("SampleMetadata" = rownames(dlist[[2]]), dlist[[2]]), file = paste0(Results.path, pref, "SampleMetadata.csv"), sep = "\t", quote = F, row.names = F)
  write.table(data.frame("VariableMetadata" = rownames(dlist[[3]]), dlist[[3]]), file = paste0(Results.path, pref, "VariableMetadata.csv"), sep = "\t", quote = F, row.names = F)
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
find.limits <- function(x,
                        y) {
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
#' @param palpha Points transparency (between 0 and 1)
#' @param psize Points size
#' @keywords scales, ggplot
#' @return ggplot list of aestethic
#' @export
#' @examples
#' plotheme.auto()
plotheme.auto <- function(Samples.grp = NULL,
                          Variables.grp = NULL,
                          limits = NULL,
                          Legend.L = T,
                          colorL = F,
                          labels = list(title = "", x = "", y = ""),
                          geom_path = F,
                          labelsL = F,
                          palpha = 0.8,
                          psize = 0.8) {
  require(ggplot2)
  opls.ggplotheme.auto <- list(
    if(!is.null(Samples.grp) & geom_path == T){geom_path(alpha = 0.4)},
    if(!is.null(Samples.grp)){labs(colour = Samples.grp)},
    if(!is.null(Samples.grp)){aes(group = as.factor(get(Samples.grp)), color = as.factor(get(Samples.grp)))},
    if(!is.null(Variables.grp)){aes(color = as.factor(get(Variables.grp)))},
    if(!is.null(Variables.grp)){labs(color = Variables.grp)},
    if(colorL == F) {scale_colour_grey()},
    geom_hline(yintercept = 0, linetype = 2, color = "grey"),
    geom_vline(xintercept = 0, linetype = 2, color = "grey"),
    geom_point(alpha = palpha, size = psize),
    if(!is.null(limits)){xlim(limits[1,1], limits[1,2])},
    if(!is.null(limits)){ylim(limits[2,1], limits[2,2])},
    ggplot_SD.theme,
    labs(labels),
    if(labelsL == T){geom_text(vjust = -0.8)},
    if(Legend.L == F){theme(legend.position = 0)} else {theme(legend.position = c(0.005,0.005), legend.justification = c(0,0), legend.direction = "horizontal", legend.title = element_blank())}
  )
  return(opls.ggplotheme.auto)
}


#' Overview of 3 level datas
#'
#' Perform a quick analysis of a 3 levels list by showing for each variables : the average value,
#' the Coefficient of Variation, the skew and the percentage of zero values. Calculations is done by
#' grouping factors for each variables and returned in a list with the plot.
#' @param var2.names vector of grouping factors names
#' @param IC Logical for IC95 calculation (could be slow if true)
#' @param plotL Logical for generating plot
#' @param alpha Transparency value of graphs
#' @param labels list of graph labels (title, x, y)
#' @param x.labL Logical for showing point labels
#' @inheritParams dlist.subset
#' @keywords summary, graph
#' @return Return a list with the calculations results as a data frame and the plot if asked.
#' @export
#' @examples
#' dlist.summary()
dlist.summary <- function(dlist,
                               var2.names = NULL,
                               IC = F,
                               plotL = F,
                               alpha = 0.2,
                               labels = list(title = "", x = "", y = ""),
                               x.labL = F) {
  require(reshape2) ; require(plyr) ; require(dplyr)
  check.list.format(dlist)
  temp.data.0 <- merge(dlist[[2]][var2.names], dlist[[1]], by = 0)
  temp.data <- melt(temp.data.0[-1], id.vars = var2.names, value.name = "value", variable.name = "variable")
  if(IC == T){
    temp.data.summary <- ddply(temp.data, c(var2.names, "variable"), summarise,
                               "N" = length(value),
                               "NA" = length(which(is.na(temp.data[[1]]) == T)),
                               "Zero" = length(which(value == 0)),
                               "Perc_Zero" = round(Zero * 100 / N, 3),
                               "Avg" = round(mean(value, na.rm = T), 3),
                               "Median" = round(median(value, na.rm = T), 3),
                               "Sum" = round(sum(value, na.rm = T), 3),
                               "Min" = round(min(value, na.rm = T), 3),
                               "Max" = round(max(value, na.rm = T), 3),
                               "SD" = round(sd(value, na.rm = T), 3),
                               "CV" = round(SD*100/Avg, 3),
                               "IC95_min_manual" = round(Avg-2*(SD/N), 3),
                               "IC95_max_manual" = round(Avg+2*(SD/N), 3),
                               "Skew" = round(e1071::skewness(value, na.rm = T), 3),
                               "Avg_IC95_min" = round(t.test(na.omit(value))$conf.int[1], 3),
                               "Avg_IC95_max" = round(t.test(na.omit(value))$conf.int[2], 3),
                               "Quant1" = round(quantile(value, na.rm = T)[1], 3),
                               "Quant2" = round(quantile(value, na.rm = T)[2], 3),
                               "Quant3" = round(quantile(value, na.rm = T)[3], 3),
                               "Quant4" = round(quantile(value, na.rm = T)[4], 3)
    )
  } else {
    temp.data.summary <- ddply(temp.data, c(var2.names, "variable"), summarise,
                               "N" = length(value),
                               "NA" = length(which(is.na(temp.data[[1]]) == T)),
                               "Zero" = length(which(value == 0)),
                               "Perc_Zero" = round(Zero * 100 / N, 3),
                               "Avg" = round(mean(value, na.rm = T), 3),
                               "Median" = round(median(value, na.rm = T), 3),
                               "Sum" = round(sum(value, na.rm = T), 3),
                               "Min" = round(min(value, na.rm = T), 3),
                               "Max" = round(max(value, na.rm = T), 3),
                               "SD" = round(sd(value, na.rm = T), 3),
                               "CV" = round(SD*100/Avg, 3),
                               "IC95_min_manual" = round(Avg-2*(SD/N), 3),
                               "IC95_max_manual" = round(Avg+2*(SD/N), 3),
                               "Skew" = round(e1071::skewness(value, na.rm = T), 3)
    )
  }
  if(plotL == T) {
    require(ggplot2)
    # Data prep
    temp.data.1 <- subset(temp.data.summary, select = c("variable", "Avg", "CV", "Skew", "Perc_Zero"))
    temp.plot <- melt(temp.data.1, variable.name = "Measure")
    temp.plot$variable <- factor(temp.plot$variable, levels = levels(reorder(temp.data.summary$variable, temp.data.summary$Avg)))
    # plot variables
    x <- "variable"
    y <- "value"
    group <- "Measure"
    yline <- rbind(data.frame("Measure" = "Skew", "yint" = c(0,-1,1), "color_var" = c("black", "red", "blue")),
                   data.frame("Measure" = "Perc_Zero", "yint" = c(50, 0, 100), "color_var" = c("black", "red", "blue")))
    # plot
    plot1 <- ggplot(temp.plot, aes(x = get(x), y = get(y), ymin = 0, ymax = get(y), color = get(group))) +
      geom_hline(data = yline, aes(yintercept = yint), color = yline$color_var, linetype = 2, alpha = alpha) +
      geom_pointrange(alpha = alpha, size = 0.02) +
      facet_grid(paste0(group, "~."), scale = "free_y") +
      ggplot_SD.theme +
      theme(legend.position = 0) +
      labs(labels)
    if(x.labL == F){
      plot1 <- plot1 + ggplot_SD_nox_lab
    } else {
      plot1 <- plot1 + ggplot_SD_lab90
    }
    return(list("Data" = temp.data.summary, "Plot" = plot1))
  }
  return(list("Data" = temp.data.summary))
}


#' dlist.summary.2
#'
#' Description of the function
#' @param var2names List of grouping factor, if blank calculation while be done for eache variables.
#' @inheritParamsdlist.summary
#' @inheritParams dlist.subset
#' @keywords x1, x2, x3
#' @return result of the function
#' @export
#' @examples
#' dlist.summary.2()
dlist.summary.2 <- function(dlist,
                            var2names = NULL,
                            plotL = F,
                            alpha = 0.2,
                            labels = list(title = "", x = "", y = ""),
                            x.labL = F){
  require(dplyr) ; require(tidyr) ; require(e1071) ## Load packages
  check.list.format(dlist)
  dlist <- lapply(dlist, tbl_df) ## transform data to tbl class
  if(length(which(var2names %in% c("N", "NA", "Zero", "Perc_Zero", "Avg", "Median",
                                   "Sum", "Min", "Max", "SD", "CV", "IC95_min_manual",
                                   "IC95_max_manual", "Skew") == T)) > 0) {
    stop("var2names conflict with calculated variable. Must be different than : N, NA, Zero, Perc_Zero, Avg, Median, Sum, Min, Max, SD, CV, IC95_min_manual, IC95_max_manual, Skew")}
  temp.summary <- bind_cols(dlist[[2]], dlist[[1]]) %>%
    gather_("variable", "value", names(dlist[[1]])) %>%
    group_by_(.dots = c(var2names, "variable")) %>%
    summarise("N" = length(value),
              "NA" = length(which(is.na(value[[1]]) == T)),
              "Zero" = length(which(value == 0)),
              "Perc_Zero" = round(Zero * 100 / N, 3),
              "Avg" = round(mean(value, na.rm = T), 3),
              "Median" = round(median(value, na.rm = T), 3),
              "Sum" = round(sum(value, na.rm = T), 3),
              "Min" = round(min(value, na.rm = T), 3),
              "Max" = round(max(value, na.rm = T), 3),
              "SD" = round(sd(value, na.rm = T), 3),
              "CV" = round(SD*100/Avg, 3),
              "IC95_min_manual" = round(Avg-2*(SD/N), 3),
              "IC95_max_manual" = round(Avg+2*(SD/N), 3),
              "Skew" = round(e1071::skewness(value, na.rm = T), 3))
  if(plotL == T) {
    require(ggplot2)
    temp.plot <- select(ungroup(temp.summary), variable, Avg, CV, Skew, Perc_Zero) %>%
      arrange(Avg) %>%
      mutate(variable = factor(variable, variable)) %>%
      gather("Measure", "value", -variable) %>%
      group_by(Measure)
    yline <- rbind(data.frame("Measure" = "Skew", "yint" = c(0,-1,1), "color_var" = c("black", "red", "blue")),
                   data.frame("Measure" = "Perc_Zero", "yint" = c(50, 0, 100), "color_var" = c("black", "red", "blue")))
    plot1 <- ggplot(temp.plot, aes(variable, value, ymin = 0, ymax = value, color = Measure)) +
      geom_hline(data = yline, aes(yintercept = yint), color = yline$color_var, linetype = 2, alpha = alpha) +
      geom_pointrange(alpha = alpha, size = 0.02) +
      facet_grid(Measure~., scale = "free_y") +
      ggplot_SD.theme +
      theme(legend.position = 0) +
      labs(labels)
    if(x.labL == F){
      plot1 <- plot1 + ggplot_SD_nox_lab
    } else {
      plot1 <- plot1 + ggplot_SD_lab90
    }
    return(list("Data" = temp.summary, "Plot" = plot1))
  }
  return(list("Data" = temp.summary))
}


#' Perform PCA and custom plot
#'
#' Perform a PCA analysis on a 3 levels list and create custom plot.
#' @param Samples.grp Factor name used for grouping samples (affect plot only)
#' @param Variables.grp Factor name used for grouping variables (affect plot only)
#' @param Legend.L Logical for drawing legends
#' @param colorL Logical for using color or greyscale
#' @param Samp.lab.L Logical for drawing Sample labels
#' @param Var.lab.L Logical for drawing Variable labels
#' @param palpha Transparency for loadings points (between 0 and 1)
#' @param psize Loadings points size
#' @inheritParams densplot
#' @inheritParams dlist.subset
#' @keywords pca, ggplot
#' @return A list with [1] pca results, [2] plot as grobs.
#' @export
#' @examples
#' dlist.pca()
dlist.pca <- function (dlist,
                       Samples.grp = NULL,
                       Variables.grp = NULL,
                       Legend.L = F,
                       colorL = F,
                       Samp.lab.L = F,
                       Var.lab.L = T,
                       palpha = 0.8,
                       psize = 0.8,
                       ShowPlot = T) {
  require(ropls) ; require(ggplot2) ; require(gridExtra)
  check.list.format(dlist)
  temp.pca <- opls(dlist[[1]], predI = 2, plotL = F)
  temp.scores <- merge(dlist[[2]], data.frame(temp.pca$scoreMN),
                       by.x = 0, by.y = 0)
  limits1 <- find.limits(temp.scores$p1, temp.scores$p2)
  temp.loadings <- merge(dlist[[3]], data.frame(temp.pca$loadingMN),
                         by.x = 0, by.y = 0)
  limits2 <- find.limits(temp.loadings$p1, temp.loadings$p2)
  labels1 <- list(title = paste0("Scores plot ", temp.pca$descriptionMC[1],
                                 " samples\n(", temp.pca$descriptionMC[4], " missing values)"),
                  x = paste0("p1 (", temp.pca$modelDF$R2X[1] * 100, " %)"),
                  y = paste0("p2 (", temp.pca$modelDF$R2X[2] * 100, " %)"))
  labels2 <- list(title = paste0("Loadings plot\n", temp.pca$descriptionMC[2],
                                 " variables (", temp.pca$descriptionMC[3], " excluded)"),
                  x = paste0("p1 (", temp.pca$modelDF$R2X[1] * 100, " %)"),
                  y = paste0("p2 (", temp.pca$modelDF$R2X[2] * 100, " %)"))
  plot1 <- ggplot(temp.scores, aes(p1, p2)) + plotheme.auto(Samples.grp,
                                                            Variables.grp = NULL, limits = limits1, Legend.L, colorL,
                                                            labels1, geom_path = T, labelsL = Samp.lab.L)
  plot2 <- ggplot(temp.loadings, aes(p1, p2, label = Row.names)) +
    plotheme.auto(Samples.grp = NULL, Variables.grp, limits = limits2,
                  Legend.L, colorL, labels1, geom_path = F, labelsL = Var.lab.L, palpha = palpha, psize = psize)
  if(ShowPlot == T) {
    return(list(PCA = temp.pca, Plot = grid.arrange(plot1, plot2, nrow = 1)))
  } else {
    return(list(PCA = temp.pca, Plot = arrangeGrob(plot1, plot2, nrow = 1)))
  }
}


#' Transform a datamatrix and replace zero
#'
#' Do any of the two following function : Replace zero by the minimum value divided by 2 and/or tranform data
#' with Trans.fun function (by default log2).
#' @param data a dataframe with value only to transform
#' @param TransL Should transformation be made
#' @param Trans.fun Which transformation to perform (default log2)
#' @param ZeroL Should zero be replaced by minimum value divided by 2
#' @param ... Argument to pass to Trans.fun
#' @keywords transform
#' @return The resulting datamatrix only
#' @export
#' @examples
#' dataframe.transform()
dataframe.transform <- function(data,
                                ZeroL = T,
                                TransL = F,
                                Trans.fun = log2,
                                ...) {
  if(isTRUE(ZeroL)) {temp.data <- data.frame(apply(data, 2, function(x) {x[x == 0] <- min(x[x!=0], na.rm = T)/2 ; x}))}
  if(isTRUE(TransL)) {return(Trans.fun(data, ...))}
  return(temp.data)
}


#' Transform function for dlist
#'
#' Replace zero by the minimum value / 2 and if checked, tranform data with Trans.fun function (by default log2).
#' Same function as dataframe.transform(), but formatted for direct use on dlist.
#' @param ... Argument passed to dataframe.transform()
#' @inheritParams dlist.subset
#' @keywords transform, dlist
#' @return The resulting datamatrix only
#' @export
#' @examples
#' dlist.transform()
dlist.transform <- function(dlist,
                            ...) {
  check.list.format(dlist)
  dlist[[1]] <- dataframe.transform(dlist[[1]], ...)
  return(dlist)
}


#' Perform OPLS and custom plots
#'
#' Perform an OPLS analysis using "ropls" package and plot results if asked (scores, loadings, VIPs & summary).
#' @param Opls.y Name of the grouping factor used in OPLS
#' @param Samples.grp Vector of samples grouping factor (for plot color only)
#' @param Variables.grp Vector of variable grouping factor (for plot color only)
#' @param Legend.L Logical to draw legends
#' @param colorL Logical to use colors instead of black & white
#' @param LabelsL Logical to show labels on plot
#' @param ShowPlot
#' @inheritParams densplot
#' @inheritParams dlist.subset
#' @keywords opls, ggplot
#' @return A list with (OPLS, Plot, VIPS)
#' @export
#' @examples
#' dlist.opls()
dlist.opls <- function (dlist,
                        Opls.y,
                        Samples.grp = NULL,
                        Variables.grp = NULL,
                        Legend.L = T,
                        colorL = F,
                        VIP.thr = 1,
                        LabelsL = F,
                        ShowPlot = T) {
  require(ropls) ; require(ggplot2) ; require(gridExtra)
  check.list.format(dlist)
  temp.opls <- opls(dlist[[1]], dlist[[2]][[Opls.y]], orthoI = NA, plotL = F)
  temp.scores <- merge(dlist[[2]], data.frame(temp.opls$scoreMN), by.x = 0, by.y = 0)
  temp.scores <- merge(temp.scores, data.frame(temp.opls$orthoScoreMN), by.x = "Row.names", by.y = 0, all = T)
  limits1 <- find.limits(temp.scores$p1, temp.scores$o1)
  temp.loadings <- merge(dlist[[3]], data.frame(temp.opls$loadingMN), by.x = 0, by.y = 0)
  temp.loadings <- merge(temp.loadings, data.frame(temp.opls$orthoLoadingMN), by.x = "Row.names", by.y = 0, all = T)
  limits2 <- find.limits(temp.loadings$p1, temp.loadings$o1)
  ## var
  labels1 <- list(title = paste0("Scores plot ", temp.opls$descriptionMC[1], " samples\n(", temp.opls$descriptionMC[4], " missing values)"),
                  x = paste0("pred. comp. 1 of ", Opls.y, " (", temp.opls$modelDF$R2X[1]*100, " %)"),
                  y = paste0("o1 (", temp.opls$modelDF$R2X[2]*100, " %)"))
  labels2 <- list(title = paste0("Loadings plot\n", temp.opls$descriptionMC[2], " variables (", temp.opls$descriptionMC[3], " excluded)"),
                  x = paste0("p1 (", temp.opls$modelDF$R2X[1]*100, " %)"),
                  y = paste0("o1 (", temp.opls$modelDF$R2X[2]*100, " %)"))
  ## plots
  plot1 <- ggplot(temp.scores, aes(p1, o1)) +
    plotheme.auto(Samples.grp, Variables.grp = NULL, limits = limits1, Legend.L, colorL, labels1, geom_path = T, labelsL = LabelsL)
  plot2 <- ggplot(temp.loadings, aes(p1, o1, label = Row.names)) +
    plotheme.auto(Samples.grp = NULL, Variables.grp, limits = limits2, Legend.L, colorL, labels1, geom_path = F, labelsL = LabelsL)

  ##
  temp.VIPs <- ggplot_opls(temp.opls, VIP.thr = VIP.thr, ShowPlot = F)
  ## Result
  if(ShowPlot == T) {
    return(list("OPLS" = temp.opls,
                "Plots" = list("ScoresPlot" = plot1,
                               "LoadingsPlot" = plot2,
                               "VipsPlot" = temp.VIPs$Plot,
                               "SummaryTable" = tableGrob(temp.opls$summaryDF, theme = ttheme_minimal(base_size = 10, padding = unit(c(2,2), "mm")))),
                "DrawPlots" = grid.arrange(plot1,
                                           plot2,
                                           tableGrob(temp.opls$summaryDF, theme = ttheme_minimal(base_size = 10, padding = unit(c(2,2), "mm"))),
                                           temp.VIPs$Plot,
                                           nrow = 2, heights = c(3,1)),
                "VIPs" = temp.VIPs$VIPs))
  } else {
    return(list("OPLS" = temp.opls,
                "Plots" = list("ScoresPlot" = plot1,
                               "LoadingsPlot" = plot2,
                               "VipsPlot" = temp.VIPs$Plot,
                               "SummaryTable" = tableGrob(temp.opls$summaryDF, theme = ttheme_minimal(base_size = 10, padding = unit(c(2,2), "mm")))),
                "DrawPlots" = arrangeGrob(plot1,
                                          plot2,
                                          tableGrob(temp.opls$summaryDF, theme = ttheme_minimal(base_size = 10, padding = unit(c(2,2), "mm"))),
                                          temp.VIPs$Plot,
                                          nrow = 2, heights = c(3,1)),
                "VIPs" = temp.VIPs$VIPs))
  }
}


#' 2 dimension plot with density
#'
#' Draw a 2 dimensions plot with density graph on each axis
#' @param Data Long format data to plot
#' @param x Column name for x values as string
#' @param y Column name for x values as string
#' @param group Column name for grouping factor as string
#' @param color Point color (default is "black")
#' @param alpha alpha value of point (default is 0.5)
#' @param labels List for labels with (title, x, y)
#' @param ShowPlot Logical to directly draw the plot or return the grobs (faster)
#' @keywords ggplot
#' @return a ggplot
#' @export
#' @examples
#' densplot()
densplot <- function(Data,
                     x,
                     y,
                     group = NULL,
                     color = "black",
                     alpha = 0.5,
                     labels = list(title = "", x = "", y = ""),
                     ShowPlot = T) {
  plot1 <- ggplot(Data, aes_string(x = x, fill = group)) +
    geom_density(adjust = 1/5, alpha = 0.25, color = color) +
    labs(list(title = labels[[1]], x = "", y = "density")) +
    theme_bw()
  plot2 <- ggplot(Data, aes_string(x = y, fill = group)) +
    geom_density(adjust = 1/10, alpha = 0.25, color = color) +
    labs(list(title = "", x = "", y = "density")) +
    theme_bw() +
    coord_flip()
  plot3 <- ggplot(Data, aes_string(x = x, y = y)) +
    geom_point(size = 0.5, aes_string(fill = group), alpha = alpha) +
    labs(list(title = "", x = labels[[2]], y = labels[[3]])) +
    theme_bw()
  if(any(!is.null(group) & group %in% names(Data)) == T) {
    plot4 <- g_legend(plot1)
  } else {
    plot4 <- ggplot() + theme(plot.background = element_blank(), panel.background = element_blank()) + theme(legend.position = "none")
  }
  plot <- arrangeGrob(plot1 + theme(legend.position = "none"),
               plot4,
               plot3 + theme(legend.position = "none"),
               plot2 + theme(legend.position = "none"),
               ncol = 2,
               heights = c(0.7,2),
               widths = c(2,0.7))
  if(ShowPlot == T) {return(plot(plot))} else {return(plot)}
}


#' 3 levels t-test
#'
#' Perform t-test on a 3 levels list between groups
#' @param group The grouping factor
#' @inheritParams dlist.subset
#' @keywords t-test
#' @return A dataframe with p.valus, x and y estimates, delta mean and CV
#' @export
#' @examples
#' data.ttest()
data.ttest <- function(dlist,
                       group) {
  check.list.format(dlist)
  ## Data prep
  temp.data <- melt(as.matrix(dlist), varnames = c("sample", "variable"))
  temp.data <- merge(temp.data, group, by.x = "sample", by.y = 0)
  print("Data OK, performing t-test")
  Statistic.summary <- ddply(temp.data, c("variable"), function(x) {
    test.try <- try(t.test(as.formula(paste("value ~",names(group))), data = x))
    if (class(test.try) != "try-error") {
      w <- t.test(as.formula(paste("value ~", names(group))), data = x)
      with(w, data.frame(statistic,
                         p.value,
                         "x.estimate" = estimate[1],
                         "y.estimate" = estimate[2],
                         "mean.delta" = abs(estimate[1] - estimate[2]),
                         "CV.mean.delta" = abs(estimate[1] - estimate[2])/max(c(estimate[1], estimate[2]))))
    }
  })
  return(Statistic.summary)
}


#' 3 Levels list from xcmsSet Object
#'
#' Create a 3 levels list from the result of xcms integration.
#' @param x xcmsSet Object
#' @keywords 3levels.list, xcmsSet
#' @return A dlist : 3 levels list (Datamatrix, SampleMetadata, VariableMetadata)
#' @export
#' @examples
#' xcmsSet.Result.List()
xcmsSet.Result.List <- function(x) {
  if(class(x) != "xcmsSet"){stop("x must be an xcmsSet object")}
  require(xcms)
  PkTable <- peakTable(x)
  Datamatrix <- data.frame(t(subset(PkTable, select = rownames(x@phenoData))))
  rownames(PkTable) <- colnames(Datamatrix)
  Variable.Metadata <- data.frame(subset(PkTable, select = !colnames(PkTable) %in% rownames(Datamatrix)))

  list("Datamatrix" = Datamatrix,
       "SampleMetadata" = x@phenoData,
       "VariableMetadata" = Variable.Metadata)
}




#' Plot OPLS results with VIPs
#'
#' Description of the function
#' @param data Result of ropls with opls method
#' @param VIP.thr VIP threshold for groups (1 by default)
#' @inheritParams densplot
#' @keywords ggplot, opls
#' @return Return a list with Plot and VIPs
#' @export
#' @examples
#' ggplot_opls()
ggplot_opls <- function(data, VIP.thr = 1, ShowPlot = T) {
  # Test
  if(class(data) != "opls"){print("Optimized for ropls::opls resulting object") ; stop()}
  if(!data$typeC %in% c("OPLS", "OPLS-DA")){print("Optimized for OPLS-DA results") ; stop()}
  # Data prep
  temp.plot <- data.frame("Vip" = data$vipVn)
  pos.list  <- subset(data.frame(data$loadingMN), p1 > 0)
  neg.list  <- subset(data.frame(data$loadingMN), p1 < 0)
  temp.plot[rownames(temp.plot) %in% rownames(neg.list),] <- temp.plot[rownames(temp.plot) %in% rownames(neg.list),] * -1
  ## Plot var
  plot.data <- temp.plot
  x <- rownames(temp.plot)
  y <- "Vip"
  VIP.subset<- subset(temp.plot, abs(Vip) >= VIP.thr)
  Select.var.VIP <- rownames(VIP.subset)
  labels <- list(title = "", x = paste0("Variables (", data$descriptionMC[2], " of which ", length(Select.var.VIP), " have a VIP > ", VIP.thr, ")"), y = "Score VIP")
  # Plot
  plot1 <- ggplot(plot.data, aes(x = reorder(x, get(y)), y = get(y), ymin = 0, ymax = get(y))) +
    geom_hline(yintercept = c(-1,-2,-VIP.thr, VIP.thr, 2), alpha = 0.4, linetype = 2) +
    geom_pointrange(alpha = 0.6, size = 0.1, aes(color = abs(get(y)) >= VIP.thr)) +
    labs(labels) +
    ggplot_theme_sly +
    ggplot_SD_nox_lab +
    theme(legend.position = 0)
  if(ShowPlot == T) {
    return(list(VIPs = VIP.subset, Plot = grid.arrange(plot1)))
  } else {
    return(list(VIPs = VIP.subset, Plot = arrangeGrob(plot1)))
  }
}


#' TITLE
#'
#' Description of the function
#' @param x Argument description
#' @keywords x1, x2, x3
#' @return result of the function
#' @export
#' @examples
#' test()
test <- function(x) {
  print(x)
}





require(ggplot2)
ggplot_SD_lab90 <- theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust = 0.95))
ggplot_labels_strip_0 <-   theme(strip.text.y = element_text(angle=0))
ggplot_SD_nox_lab <- theme(axis.text.x = element_blank(), axis.ticks = element_blank())
ggplot_SD_noy_lab <- theme(axis.text.y = element_blank(), axis.ticks = element_blank())
ggplot_no_labels <- ggplot_SD_nox_lab + ggplot_SD_noy_lab
ggplot_theme_sly <- theme_classic() +  theme(axis.line.x = element_line(), axis.line.y = element_line())
ggplot_SD.theme <- theme_bw() + theme(panel.grid.minor = element_blank(),
                                      panel.grid.major = element_blank())

