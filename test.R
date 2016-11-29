## DATA
require(tidyverse) ; require(data.table) ; require("SDjoygret") ; require(tibble) ; require(dplyr)
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
lapply(data[[3]], tbl_df)
lapply(data[[2]], data.table)


require(dtplyr)

temp <- rbenchmark::benchmark(
  dplyr.tibble      = filter(data[[2]][[1]], RowID %in% c("HU_015", "HU_011")),
  dplyr.d.t         = filter(data[[3]][[1]], RowID %in% c("HU_015", "HU_011")),
  data.table.d.t    = data[[3]][[1]][RowID %in% c("HU_015", "HU_011")],
  replications = 1000
)







### DATA.TABLE TUTORIAL
# https://cran.r-project.org/web/packages/data.table/vignettes/datatable-intro.html
flights <- fread("https://raw.githubusercontent.com/wiki/arunsrinivasan/flights/NYCflights14/flights14.csv")

## Select and order rows
flights[1:10]
flights[month == 7 & origin == "JFK"]
flights[order(month, -flight)][1:10]

## Select column
flights[,.(arr_delay)] ## column
flights[,arr_delay]    ## vector
flights[,.(month, day, cancelled)] ## multiples
flights[,.(delay_arr = arr_delay, delay_dep = dep_delay)] ## subset and rename columns

## calcul on column
flights[,sum((arr_delay + dep_delay) > 0)]
flights[origin == "JFK" & month == 7, length(dest)] ## count row
flights[origin == "JFK" & month == 7, .N]           ## count row (any column) : .N holds the number of observations in the current group

## calcul
flights[origin == "JFK" & month == 6,.(Av_dep_delay = mean(dep_delay), Av_ar_delay = mean(arr_delay))] ## on rows subset

## use names as reference with arguments : with = F
flights[, c("arr_delay", "dep_delay"), with = F]
flights[, !c("arr_delay", "dep_delay"), with = F]
flights[, -c("arr_delay", "dep_delay"), with = F]

## by group
flights[,.(Av_dep_delay = mean(dep_delay), Av_ar_delay = mean(arr_delay)), .(origin, month)] ## by group
flights[, .N, origin] ## number of observation by origin
flights[, .N, "origin"] ## by accept quoted variables
flights[, .(.N), "origin"] ## same
flights[carrier == "AA", .(.N), .(origin, dest)] ## subset and group by
flights[carrier == "AA", .(Av_arr_delay = mean(arr_delay), Av_dep_delay = mean(dep_delay)), .(origin, dest, month)] ## subset and group by
flights[carrier == "AA", .(Av_arr_delay = mean(arr_delay), Av_dep_delay = mean(dep_delay)), .(origin, dest, month)][order(origin, dest, month)] ## subset and group by with order, but better use keys :
flights[carrier == "AA", .(Av_arr_delay = mean(arr_delay), Av_dep_delay = mean(dep_delay)), keyby = .(origin, dest, month)] ## keys are a lot faster

## chaining
flights[carrier == "AA", .N, by = .(origin, dest)][order(origin, -dest)] ## horizontally
flights[carrier == "AA", .N, by = .(origin, dest)                        ## or vertically
        ][order(origin, -dest)]

## group by functions
flights[, .N, .(dep_delay>0, arr_delay>0)]

## group by and return mean of all column (summarize_all)
flights[carrier == "AA",                         ## select one carrier
        lapply(.SD, mean),                       ## compute mean of all column (.SD)
        .(origin, dest, month),                  ## group by origin, dest and month
        .SDcols = c("arr_delay", "dep_delay")]   ## only include this two column in .SD (accept inverse : all but... (- or !)

# the same :
flights[carrier == "AA",                         ## select one carrier
        lapply(.(arr_delay, dep_delay), mean),   ## compute mean of all column (.SD)
        .(origin, dest, month)]                  ## only include this two column in .SD

## concatenate two column
flights[,.(groups = paste0(carrier, flight), carrier, flight), by = origin]

## concatenate two column : all flights of a carrier by origin
flights[,.(Flights_names = paste(c(origin, flights))), carrier]




