% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Package.1_functions.R
\name{import.xls.2}
\alias{import.xls.2}
\title{Import Excel's sheets to a list}
\usage{
import.xls.2(Data.path)
}
\arguments{
\item{Data.path}{Path to the excel file to import}
}
\value{
A list of dataframe
}
\description{
Import each sheet of an excel file to a list as a data.frame. This function use the readxl package, which seems
to handle large files.
}
\examples{
import.xls.2()
}
\keyword{excel}
\keyword{import,}
\keyword{list,}
