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

sacurine.dt <- data[[3]]
save(sacurine.dt, file = "./data/sacurine.dt.rda")

## DF with rownames
dataframe <- data.frame(A = 1:10,
           B = 1:10,
           C = 1:10,
           row.names = paste0("A", 1:10))

## check
lapply(data, check.list.format)
dlist <- data[[1]]


## DEVELOPMENT



















