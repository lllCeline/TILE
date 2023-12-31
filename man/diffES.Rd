% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diffES.R
\name{diffES}
\alias{diffES}
\title{limma analysis of cell type score}
\usage{
diffES(score_matrix, groupinfo)
}
\arguments{
\item{score_matrix}{a data matrix with rows labeling features immune cell type, and columns labeling samples.}

\item{groupinfo}{a matrix with two columns: 'sample' column labeling sample id,'groups' column labeling group information}
}
\description{
Different immune cell types were analyzed by the limma package
}
\examples{
result <- diffES(score_matrix = cell_enrich,groupinfo = demo_groupinfo)
}
