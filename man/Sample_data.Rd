\name{Sample_data}
\alias{Sample_data}
\docType{data}
\title{
Sample data for analysis
}
\description{
Sample datasets for analysis using HuberNet and Net(TF, PWG, Annotation); Sample datasets for analysis using HuberLasso, Lasso, HuberENET, and ENET (TF, PWG).
}
\usage{APGD::Sample_data}
\format {Datasets with 3 dataframes, which are:
}
\arguments{

\item{TF} {Expression level for one TF with 50 sample size.}
\item{PWG} {Expression levels for 131 PWGs with 50 sample size.}
\item{Annotation} {Annotation file for 131 PWGs, which are included in 5 Pathways. (Only for HuberNet and Net methods)}

}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
Sample_data <- APGD::Sample_data

X <- Sample_data$PWG

y <- Sample_data$TF

Annotation <- Sample_data$Annotation
}
\keyword{datasets}
