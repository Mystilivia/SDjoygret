## DATA
require(tidyverse) ; require(data.table)
data <- list("data.frame" = list("Datamatrix"       = data.frame("RowID" = rownames(sacurine[[1]]), sacurine[[1]]),
                                 "SampleMetadata"   = data.frame("RowID" = rownames(sacurine[[2]]), sacurine[[2]]),
                                 "VariableMetadata" = data.frame("RowID" = rownames(sacurine[[3]]), sacurine[[3]])),
             "tibble"     = list("Datamatrix"       = tbl_df(data.frame("RowID" = rownames(sacurine[[1]]), sacurine[[1]])),
                             "SampleMetadata"       = tbl_df(data.frame("RowID" = rownames(sacurine[[2]]), sacurine[[2]])),
                             "VariableMetadata"     = tbl_df(data.frame("RowID" = rownames(sacurine[[3]]), sacurine[[3]]))),
             "data.table" = list("Datamatrix"       = data.table("RowID" = rownames(sacurine[[1]]), sacurine[[1]]),
                                 "SampleMetadata"   = data.table("RowID" = rownames(sacurine[[2]]), sacurine[[2]]),
                                 "VariableMetadata" = data.table("RowID" = rownames(sacurine[[3]]), sacurine[[3]]))
             )

## check
lapply(data, dlist.class)



check.list.format <- function (dlist, to.data.table.L = F, return.dlist = T) {
  require("data.table")
  require("tidyverse")
  if(!is.list(dlist)){stop("Data should be a list with (1) Datamatrix (2) Sample.Metadata (3) Variable.Metadata")}
  temp.data.str <- dlist.class(dlist)
  if(any(temp.data.str[, class.d.t] == F) & any(temp.data.str[,class.d.f] == F)) {stop("List levels should be data.frame or data.table")}
  if(!any(temp.data.str[, class.d.t] == F)) { ## all are data.table
    if(!identical(names(dlist[[1]])[-1], dlist[[3]][[1]])){stop("Datamatrix colnames should be identical of Variable.Metadata rownames")}
    if(!identical(dlist[[1]][[1]], dlist[[2]][[1]])){stop("Datamatrix rownames should be identical of Sample.Metadata rownames")}
    dim.temp <- lapply(dlist, dim)
    if(!dim.temp[[1]][1] == dim.temp[[2]][1]){stop("Datamatrix row number should be the same as Sample.Metadata")}
    if(!dim.temp[[1]][2]-1 == dim.temp[[3]][1]){stop("Datamatrix col number should be the same as Variable.Metadata row number")}
    print("Data format is ok")
    print(temp.data.str)
    if(return.dlist == T) {return(dlist)}
  } else {
    if(!any(temp.data.str[,class.d.f] == F)) { ## all are data.table
      if(!identical(colnames(dlist[[1]]), rownames(dlist[[3]]))){stop("Datamatrix colnames should be identical of Variable.Metadata rownames")}
      if(!identical(rownames(dlist[[1]]), rownames(dlist[[2]]))){stop("Datamatrix rownames should be identical of Sample.Metadata rownames")}
      dim.temp <- lapply(dlist, dim)
      if(!dim.temp[[1]][1] == dim.temp[[2]][1]){stop("Datamatrix row number should be the same as Sample.Metadata")}
      if(!dim.temp[[1]][2] == dim.temp[[3]][1]){stop("Datamatrix col number should be the same as Variable.Metadata row number")}
      print("Data seems OK but are stored in data.frame")
      if(to.data.table.L == T) {
        to.data.table(dlist, T)
      } else {
        print("Data weren't converted to data.table, set arg to.data.table = T if convertion is needed")
        print(temp.data.str)
        if(return.dlist == T) {return(dlist)}
      }
    }
  }
}


require(SDjoygret)
lapply(data, dlist.class)
