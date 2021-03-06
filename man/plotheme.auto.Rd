% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Package.1_functions.R
\name{plotheme.auto}
\alias{plotheme.auto}
\title{Semi-automatic ggplot theme}
\usage{
plotheme.auto(limits = NULL, Legend.L = T, colorL = F,
  labels = list(title = "", x = "", y = ""), geom_path = F,
  geom_ellipse = F, labelsL = F, palpha = 0.8, psize = 0.8)
}
\arguments{
\item{limits}{Limits of the axis (returned by find.limits)}

\item{Legend.L}{Logical to draw legend}

\item{colorL}{Logical to use color or grey scale}

\item{labels}{list for titles (title, x, y)}

\item{geom_path}{Logical for drawing path}

\item{geom_ellipse}{Logical for drawing groups ellipses}

\item{labelsL}{Add labels to points}

\item{palpha}{Points transparency (between 0 and 1)}

\item{psize}{Points size}
}
\value{
ggplot list of aestethic
}
\description{
Returns the ploting parameters according to input, with working default for missing ones.
Used for custom ggplot of multivariate results (pca, pls, opls)
}
\examples{
plotheme.auto()
}
\keyword{ggplot}
\keyword{scales,}
