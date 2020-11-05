#' @keywords internal
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL

#' Sample amplicon sequences
#'
#' These are the three most abundant ASVs from
#' \href{https://www.authorea.com/doi/full/10.22541/au.160340221.18389016/v1}{Furneaux et al. 2020}.
#' They are derived from PacBio RSII metabarcoding of soil samples from
#' central Benin, West Africa.
#' The three sequences are from the fungal genera \emph{Penicillium}, \emph{Agaricus}
#' (included as a control), and \emph{Inocybe}.
#'
#' @format A \code{\link[Biostrings]{DNAStringSet}} with three sequences.
"tzara_sample"

#' Sample quality scores
#'
#' These are based on CCS quality scores from the first 1000 reads of long
#' amplicon PacBio RSII samples in
#' \href{https://www.authorea.com/doi/full/10.22541/au.160340221.18389016/v1}{Furneaux et al. 2020}.
#'
#' @format A \code{data.frame} with columns \code{seq_id}, \code{quality},
#' and \code{n}, giving the number of times each wuality score appears in each
#' read.
"quality_profile"
