% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Package.1_functions.R
\name{stat_table}
\alias{stat_table}
\title{Do pairwise statistical tests on dlist}
\usage{
stat_table(data, factor, value, group.by = NULL, ..., Output = c("simple",
  "complete"), debug = F)
}
\arguments{
\item{data}{A dlist}

\item{factor}{The column name with comparison}

\item{value}{The column name with values}

\item{group.by}{Grouping variable for ggpubr::compare_means}

\item{...}{Parameters to pass to ggpubr::compare_means}

\item{Ouptut}{Output mode, simple for letters list or complete for letters and p-value}
}
\value{
result of the function
}
\description{
These function take a table as entry and perform pairwise statistical test using ggpubr::compare_means for
tests and multcompView::multcompLetters to annotate test results. It returns a data.table.
}
\examples{
stat_table()
}
\keyword{dlist}
\keyword{letters,}
\keyword{significance}
\keyword{statistical}
\keyword{test,}
