% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tzara.R
\name{find_all_region_chimeras}
\alias{find_all_region_chimeras}
\title{Find chimeras in adjacent variable-conserved-variable domain triplets}
\usage{
find_all_region_chimeras(
  region_table,
  order,
  sample_column = NULL,
  read_column = "read_id",
  chimera_offset = 0,
  ...
)
}
\arguments{
\item{region_table}{(\code{data.frame} containing all the regions in
\code{order})}

\item{order}{(\code{character} with the regions to check)}

\item{sample_column}{(\code{character} scalar) extra column in region_table
which identifies different samples; if present, use
\code{\link[dada2]{isBimeraDenovoTable}} instead of
\code{\link[dada2]{isBimeraDenovo}}.}

\item{read_column}{(\code{character} scalar) column in region_table
which identifies different reads.}

\item{chimera_offset}{(\code{integer} scalar) use \code{1} if the first
region is conserved.}

\item{...}{passed to \code{find_region_chimeras}.}
}
\value{
a \code{\link[tibble]{tibble}} with columns:
    \describe{
        \item{chimset}{(\code{list} of \code{character}) the set of regions
            (usually 3) which were checked in each row.}
        \item{chims}{\code{list} of \code{character}) the sequence IDs of
            the reads which were identified as chimeric for each set of
            regions}}
}
\description{
Find chimeras in adjacent variable-conserved-variable domain triplets
}
