% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Package.1_functions.R
\name{mz.database}
\alias{mz.database}
\title{Search mz values in a database}
\usage{
mz.database(mass, database, ppm = 5)
}
\arguments{
\item{mass}{mz to search}

\item{database}{database to search into (must have 'mz' column)}

\item{ppm}{tolerance in ppm}
}
\value{
The resulting list of opls function (subsetted if min = TRUE).
}
\description{
Return all entry of database with mz between mz +- ppm
}
\examples{
mz.database(mz = 188.2500, database = data.table(ID = paste0(ID, 1:100), mz = seq(100.000, 800.000, length.out = 100)), ppm = 100)
}
\keyword{ggplot}
\keyword{opls,}
\keyword{vips,}
