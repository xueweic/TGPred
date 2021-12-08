\name{CalculateAdj}
\alias{CalculateAdj}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Calculate adjacency matrix from an annotation file.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CalculateAdj(Annotation)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Annotation}{
    n_genes*n_pathways dimensional matrix that indicates the annotation of genes within pathways information.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

  \item{Adj }{the adjacency matrix of network structure.}

}
\examples{

    Annotation <- data.matrix(Annotation)
    
    CalculateAdj(Annotation)
}



