% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tzara-package.R
\docType{package}
\name{tzara-package}
\alias{tzara}
\alias{tzara-package}
\title{tzara: Cluster long amplicons using dada2 denoising on variable regions}
\description{
To reduce computational complexity, dada2 only uses non-singletons
    as seeds for denoising.  For this strategy to work, each true sequence must
    be represented by at least two identical reads.  Especially with long
    amplicons, the probability of two reads having exactly the same errors is
    much lower than the probability of being error-free, so in practice this
    means that each true sequence must have two error-free reads.  This becomes
    problematic for rare sequences in long amplicon libraries. An alternative
    is to use hidden Markov models to cut out the most variable section of the
    targeted region and use dada2 to create denoised sequences using only that
    sequence, and then find a consensus sequence for all sequences that match
    the index region. Tzara (named after Tristan Tzara, a central figure in the
    Dada art movement) applies this method to rDNA sequences by cutting out the
    variable ITS2 region using rITSx.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://www.github.com/brendanf/tzara}
  \item Report bugs at \url{https://www.github.com/brendanf/tzara/issues}
}

}
\author{
\strong{Maintainer}: Brendan Furneaux \email{brendan.furneaux@gmail.com}

}
\keyword{internal}
