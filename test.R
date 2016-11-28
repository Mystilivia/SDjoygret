## DATA
require(tidyverse) ; require(data.table) ; require("SDjoygret")
sacurine <- SDjoygret::sacurine.dlist
data <- list("data.frame" = list("Datamatrix"       = data.frame("RowID" = sacurine[[1]][[1]], sacurine[[1]][-1]),
                                 "SampleMetadata"   = data.frame("RowID" = sacurine[[2]][[1]], sacurine[[2]][-1]),
                                 "VariableMetadata" = data.frame("RowID" = names(sacurine[[1]][-1]), sacurine[[3]][-1])),
             "tibble"     = list("Datamatrix"       = tbl_df(data.frame("RowID" = sacurine[[1]][[1]], sacurine[[1]][-1])),
                             "SampleMetadata"       = tbl_df(data.frame("RowID" = sacurine[[2]][[1]], sacurine[[2]][-1])),
                             "VariableMetadata"     = tbl_df(data.frame("RowID" = names(sacurine[[1]][-1]), sacurine[[3]][-1]))),
             "data.table" = list("Datamatrix"       = data.table("RowID" = sacurine[[1]][[1]], sacurine[[1]][-1]),
                                 "SampleMetadata"   = data.table("RowID" = sacurine[[2]][[1]], sacurine[[2]][-1]),
                                 "VariableMetadata" = data.table("RowID" = names(sacurine[[1]][-1]), sacurine[[3]][-1]))
             )

## DF with rownames
dataframe <- data.frame(A = 1:10,
           B = 1:10,
           C = 1:10,
           row.names = paste0("A", 1:10))

## check
lapply(data, check.list.format)
dlist <- data[[1]]


## DEVELOPMENT
check.list.format <- function (dlist, rownamesL = F, tibbleL = F) {
  require("tidyverse")
  if(!is.list(dlist)){stop("Data should be a list with (1) Datamatrix (2) Sample.Metadata (3) Variable.Metadata")}
  temp.data.str <- dlist.class(dlist)
  if(all(temp.data.str[,class.d.t] == F) & all(temp.data.str[,class.t] == F) & all(temp.data.str[,class.d.f] == F)) {stop("List levels should be data.frame, tibble or data.table") }
  if(rownamesL == F){
    if(!identical(dlist[[1]][[1]], dlist[[2]][[1]])){stop("Datamatrix rownames should be identical of Sample.Metadata rownames")}
    if(!identical(names(data.frame(dlist[[1]])[-1]), as.character(dlist[[3]][[1]]))){stop("Datamatrix colnames should be identical of Variable.Metadata rownames")}
  } else if(rownamesL == T){
    if(!identical(rownames(dlist[[1]]), rownames(dlist[[2]]))){stop("Datamatrix rownames should be identical of Sample.Metadata rownames")}
    if(!identical(names(dlist[[1]]), rownames(dlist[[3]]))){stop("Datamatrix colnames should be identical of Variable.Metadata rownames")}
    message("Data are well structured.")
    return(dlist)
  }
  if(all(temp.data.str[,class.t] == T)) {message("Data are stored as tibbles, well done !")
  } else if(all(temp.data.str[,class.d.t] == T)) { ## all are data.table
    if(tibbleL == F) {stop("Data are stored as data.table, to convert to tibbles : set args tibbleL to TRUE.")}
    if(tibbleL == T) {
      if(rownamesL == T) {stop("Data are stored as data.table, there shouldn't be any rownames, please check and set argument rownamesL to FALSE.")}
      if(rownamesL == F) {
        temp.list <- lapply(dlist, function(x) {tbl_df(x)})
        names(temp.list) <- names(dlist)
        message("Data were stored as data.table and were converted to tibbles.")
        return(temp.list)
      }
    }
  } else if(all(temp.data.str[,class.d.f] == T)) { ## all are data.frame
    if(tibbleL == F) {stop("Data are stored as data.frame, to convert to tibbles : set args tibbleL to TRUE.")}
    if(tibbleL == T) {
      if(rownamesL == T) {
        temp.list <- lapply(dlist, function(x) {tbl_df(rownames_to_column(x))})
        names(temp.list) <- names(dlist)
        message("Data were stored as data.frame and were converted to tibbles, rownames were added as first column.")
        return(temp.list)
      } else if(rownamesL == F) {
        temp.list <- lapply(dlist, function(x) {tbl_df(x)})
        names(temp.list) <- names(dlist)
        message("Data were stored as data.frame and were converted to tibbles.")
        return(temp.list)
      }
    }
  } else {stop("Data class isn't recognized as data.frame, tibble or data.table.")}
}




require(SDjoygret)
lapply(data, dlist.class)
