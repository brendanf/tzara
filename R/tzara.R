#' @import utils
utils::globalVariables(c(".", ".seqs", "asv.idx", "asv.seq", "border",
                         "dada.idx", "derep", "derep.idx", "derep.seq", "end",
                         "newmap", "seq.id", "start", "the_sread"))

#' @importFrom magrittr %>%
`%>%`

.onLoad <- function(libname, pkgname) {
   backports::import(pkgname)
}

#' Combine DADA2 \code{\link[dada2:derep-class]{derep}} objects into a master map
#'
#' @param dereps a (possibly named) \code{list} of (\code{\link[dada2:derep-class]{derep}} objects), or a \code{\link[tibble]{tibble}} with a column "derep" containing such a list.
#' @param .data a \code{\link[tibble]{tibble}} with the same number of rows as the length of \code{dereps}.
#' @param ... additional columns to add to the output.
#'
#' @details To be useful for further analysis, each sequence should be uniquely identified.  This can be done in several ways:
#'   \itemize{
#'     \item{\code{dereps} is a named \code{list} with unique names;}
#'     \item{\code{dereps} is a \code{\link[tibble]{tibble}} with columns (other than "derep")
#'         which uniquely identify the rows;}
#'     \item{\code{.data} is provided and its rows are unique; or}
#'     \item{\code{...} is provided and its combinations are unique.}}
#'
#' @return \code{list} with two members:
#' \describe{
#'   \item{\code{$map} (\code{\link[tibble]{tibble}})}{with columns: "\code{file}" (\code{character}), "\code{idx}" (\code{integer}), and "\code{map}" (\code{integer}), giving the mapping from the "\code{idx}"th sequence in "\code{file}" to a sequence in "\code{fasta}"}
#'   \item{\code{$fasta} (\code{\link[Biostrings]{DNAStringSet}})}{all unique sequences; the name of each sequence is an \code{integer} which matches a value in \code{map$newmap}}
#' }
#' @export
combine_derep <- function(dereps, .data = NULL, ...) {

   # handle different input types to get a tibble with the derep objects in one
   # column
   groups <- names(list(...))
   if (is.data.frame(dereps)) {
      if (length(list(...)) > 0) {
         dereps <- dplyr::bind_cols(
            tibble::tibble(.placeholder = seq_along(nrow(dereps)),
                        ...),
            dereps)
         dereps <- dplyr::select(dereps, -".placeholder")
      }
   } else {
      n <- names(dereps)
      dereps <- tibble::tibble(..., derep = dereps)
      if (!hasName(dereps, "name") && !is.null(n)) dereps$name <- n
   }

   if (!missing(.data)) {
      dereps <- dplyr::bind_cols(.data, dereps)
   }

   # Check that our derep objects are uniquely identified
   gps <- setdiff(names(dereps), "derep")
   assertthat::assert_that(dplyr::n_distinct(dereps[gps]) == nrow(dereps))

   # get all the old mappings
   # preserve the sequence names as "seq.id" if they are present
   oldmap <- dereps
   oldmap[["oldmap"]] <- lapply(oldmap[["derep"]], `[[`, "map")
   if (all(vapply(oldmap[["derep"]], assertthat::has_name, TRUE, "names"))) {
      oldmap[["seq.id"]] <- lapply(oldmap[["derep"]], `[[`, "names")
   }
   oldmap <- dplyr::select(oldmap, -derep) %>%
      tidyr::unnest() %>%
      dplyr::group_by_at(gps) %>%
      dplyr::mutate(idx = 1:dplyr::n()) %>%
      dplyr::ungroup()

   # get the old unique sequences
   olduniques <- dereps %>%
      dplyr::mutate_at("derep",
                       ~purrr::map(., .f = ~tibble::tibble(seq = names(.$uniques),
                                                           n = .$uniques))) %>%
      tidyr::unnest()

   # combine duplicate sequences among all files.
   newuniques <- olduniques %>%
      dplyr::group_by(seq) %>%
      dplyr::summarize(n = sum(n)) %>%
      dplyr::arrange(dplyr::desc(n)) %>%
      dplyr::transmute(seq = seq,
                       newmap = seq_along(seq))

   # create a mapping from the unique sequences list in each file
   # to the master unique sequence list
   newderep <- olduniques %>%
   {dplyr::left_join(
      dplyr::select(., dplyr::one_of(gps), seq) %>%
         dplyr::group_by_at(gps) %>%
         dplyr::mutate(oldmap = seq_along(seq)) %>%
         dplyr::ungroup(),
      newuniques,
      .by = "seq")}

   out <- list()
   # map from the individual sequences in each file to the master unique
   # sequence list
   out$map <- dplyr::left_join(oldmap,
                               dplyr::select(newderep, dplyr::one_of(gps),
                                             "oldmap", "newmap"),
                               by = c(gps, "oldmap")) %>%
      dplyr::select(-oldmap, map = newmap)
   #unique sequence list
   out$fasta <- Biostrings::DNAStringSet(x = purrr::set_names(newuniques$seq,
                                                              newuniques$newmap))
   class(out) <- c("multiderep", class(out))
   return(out)
}


#' Produce a map between denoised sequences and their original identifiers.
#'
#' @param derep a \code{\link[dada2]{derep-class}} object or list of such objects
#' @param dada (\link[dada2]{dada-class} object or list of such objects) the results of a call to \code{\link[dada2]{dada}} on \code{derep}
#' @param ... additional columns to add to the output.  Names included in the output by default should be avoided.
#'
#' @details Columns \code{$derep.seq} and \code{$dada.seq} contain one sequence per read
#' as plain character strings. Keeping a separate, dereplicated list of sequences and storing references to them may seem like it would be more memory efficient, but it is not necessary to do this explicitly, because this is what R already does with \code{character} vectors; only one copy of each unique string is actually kept in memory, and everything else is pointers.
#'
#' @return a \code{\link[tibble]{tibble}} with columns:
#'   \describe{
#'     \item{\code{$name} (\code{character})}{names from \code{derep}, if it is a named list of \code{\link[=add_derep_names]{named_derep}} objects. Otherwise absent unless provided in \code{...}.}
#'     \item{\code{...} (\code{character})}{any additional arguments, passed on
#'       to \code{\link[tibble]{tibble}}}
#'     \item{\code{$seq.id} (\code{character})}{the sequence identifiers from the
#'       original fasta/q file.}
#'     \item{\code{$derep.idx} (\code{integer})}{the index of each sequence in
#'       \code{derep$uniques}.}
#'     \item{\code{$derep.seq} (\code{character})}{the sequences.}
#'     \item{\code{$dada.seq} (\code{character})}{the denoised sequences.}
#' }
#' @export

dadamap <- function(derep, dada, ...) UseMethod("dadamap")
#' @export
dadamap.derep <- function(derep, dada, ...) {
   m <- tibble::tibble(seq.id = derep$names,
              derep.idx = derep$map,
              derep.seq = names(derep$uniques)[derep.idx],
              ...)

   m <- dplyr::left_join(m,
                    tibble::tibble(
                       dada.idx = dada$map,
                       derep.idx = seq_along(dada.idx),
                       dada.seq = dada$sequence[dada.idx]))
   class(m) <- c("dadamap", class(m))
   m
}

#' @export
dadamap.list <- function(derep, dada, ...) {
   assertthat::assert_that(assertthat::are_equal(length(derep), length(dada)))
   for (i in seq_along(derep)) {
      assertthat::assert_that(methods::is(derep[[i]], "derep"))
      assertthat::assert_that(methods::is(dada[[i]], "dada"))
   }

   args <- list(..., derep = derep, dada = dada)
   if (!is.null(names(derep))) {
      args[["name"]] <- names(derep)
      args <- dplyr::select(args, "name", dplyr::everything())
   }

   out <- purrr::pmap_dfr(args, dadamap.derep)
   class(out) <- c("dadamap", class(out))
   out
}


#' Individually hash biological sequences
#'
#' @param seq (\code{character} or \code{\link[Biostrings]{XStringSet}}) the
#'            sequences to hash.
#' @param algo (\code{character}) a hash algorithm supported by
#'             \code{\link[digest]{digest}}. default: "xxhash32"
#' @param len (\code{integer}) number of characters to keep from each hash
#'            string. NA (the default) to keep all characters.
#'
#' @return a \code{character} vector of the same length as \code{seq},
#'         with the hashed sequences.
#' @export
seqhash <- function(seq, algo = "xxhash32", len = NA) UseMethod("seqhash")

#' @rdname seqhash
#' @export
seqhash.character <- function(seq, algo = "xxhash32", len = NA) {
   h <- vapply(seq, digest::digest, "", algo = algo)
   if (is.na(len)) {
      return(h)
   } else {
      substring(h, 1, len)
   }
}

# @importClassesFrom dada2 derep
# setClass("named_derep", contains = "derep")

#' @rdname seqhash
#' @export
seqhash.XStringSet <- function(seq, algo = "xxhash32", len = NA) {
   seqhash.character(as.character(seq), algo = algo, len = len)
}

#' Add sequence names to a derep object
#'
#' @param derep (object of class \code{\link[dada2:derep-class]{derep}} or a \code{list} of  such objects) object(s) to add names to.
#' @param ... passed to methods
#'
#' @return (object of class \code{derep}, or a list of such objects) a
#' shallow copy of \code{derep}, with an additional member "\code{$names}",
#' giving the identifiers for the sequences from the original fasta/q file.
#' @export
add_derep_names <- function(derep, ...) UseMethod("add_derep_names")

#' @rdname add_derep_names
#' @param filename (\code{character}) name of fasta/q file that the
#'     \code{\link[dada2:derep-class]{derep}} object was derived from.
#' @export
add_derep_names.derep <- function(derep, filename, ...) {
   assertthat::assert_that(file.exists(filename))
   fqs <- ShortRead::FastqStreamer(filename, n = 1e4)
   on.exit(close(fqs))
   seqids <- character(0)
   while (length(fq <- ShortRead::yield(fqs, qualityType = "FastqQuality"))) {
      seqids <- c(seqids, as.character(fq@id))
   }
   derep[["names"]] <- seqids
   # derep <- as(derep, "named_derep")
   return(derep)
}
#' @rdname add_derep_names
#' @param filenames (\code{character}) name(s) of file(s) that a \code{list} of
#'     \code{\link[dada2:derep-class]{derep}} was derived from.
#' @export
add_derep_names.list <- function(derep, filenames = names(derep), ...) {
   if (!all(inherits(derep, "derep") |
                                  vapply(derep, is.na, TRUE))) {
      stop("length of filenames and derep do not match")
   }
   assertthat::assert_that(assertthat::are_equal(length(derep), length(names)))
   derep2 <- mapply(add_derep_names.derep, derep, filenames)
   names(derep2) <- filenames
   return(derep2)
}

#' Summarize one or more \code{\link[ShortRead:ShortReadQ-class]{ShortReadQ}} objects as a \code{\link[tibble]{tibble}}.
#'
#' @param sread (\code{\link[ShortRead:ShortReadQ-class]{ShortReadQ}} object or list of such objects) as produced by \code{\link[ShortRead]{readFastq}}.
#' @param max_ee (\code{numeric}) filter out reads with expected error greater than \code{max_ee}.  Default: Inf (no filtering)
#' @param ... (any vector) additional columns to add to the output \code{\link[tibble]{tibble}}.
#'   These should have length 1 or, if \code{sread} is a named list of \code{\link[ShortRead:ShortReadQ-class]{ShortReadQ}}, the length of \code{sread}.  Recycling behavior is as in \code{\link[tibble]{tibble}}. Avoid the name "\code{name}" if \code{sread} is a named \code{list}
#'
#' @return a \code{\link[tibble]{tibble}} with columns: \describe{
#'    \item{\code{$seq.id} (\code{character})}{sequence IDs, typically from the input fastq}
#'    \item{\code{$seq} (\code{character})}{the actual sequence reads}
#'    \item{\code{$name} (\code{character})}{the name of the source \code{\link[ShortRead:ShortReadQ-class]{ShortReadQ}}, if \code{sread} was a list. Otherwise absent.}
#'    \item{\code{...}}{any other arguments passed to summarize_sread.}}
#' If the input was empty, not a valid \code{\link[ShortRead:ShortReadQ-class]{ShortReadQ}}, or no sequences passed filtering, a \code{\link[tibble]{tibble}}
#' with zero rows is returned.
#' @export
summarize_sread <- function(sread, ..., max_ee) UseMethod("summarize_sread")

#' @export
summarize_sread.ShortReadQ <- function(sread, ..., max_ee = Inf) {
   if (!methods::is(sread, "ShortReadQ")) return(tibble::tibble(
      seq.id = character(),
      seq = character()))
   out <- tibble::tibble(
      seq.id = as.character(sread@id),
      seq = as.character(sread@sread),
      ...)
   ee <- rowSums(10^-(methods::as(sread@quality, "matrix")/10), na.rm = TRUE)
   out <- out[ee <= max_ee,,drop = FALSE]
}

#' @export
summarize_sread.list <- function(sread, ..., max_ee = Inf) {
   out = tibble::tibble(
      .seqs = lapply(sread, summarize_sread.ShortReadQ(the_sread)),
      ...)
   if (!is.null(names(sread)) && !hasName(out, "name")) {
      out$name <- names(sread)
   }
   tidyr::unnest(out, .seqs)
}

#' Combine raw and denoised reads from multiple regions of the same sequences.
#'
#' @param dadamap (\code{\link{dadamap}} object)
#' @param rawdata (\code{\link[tibble]{tibble}}) as returned by \code{\link{summarize_sread}}
#' @param key (\code{character}) no longer used??
#'
#' @details Both of the inputs should be annotated with identifying columns to
#' uniquely identify the source \code{\link[dada2:derep-class]{derep}} objects;
#' this should happen automatically if \code{\link[dada2]{derepFastq}} and
#' \code{\link[ShortRead]{readFastq}} are called on a list of filenames (there
#' will be a column "name" in the outputs of \code{\link{dadamap}} and
#' \code{\link{summarize_sread}}).
#'
#' @return a \code{\link[tibble]{tibble}} giving the sequences for each region (contents of the column given by \code{key}) for each read (uniquely identified by all other columns)
#' @export
combine_bigmaps <- function(dadamap, rawdata, key = "Region") {
      dplyr::full_join(dadamap, rawdata) %>%
      dplyr::group_by(seq.id) %>%
      dplyr::filter(any(!is.na(asv.idx))) %>%
      dplyr::mutate(seq = dplyr::coalesce(asv.seq, derep.seq, seq)) %>%
      dplyr::select(-derep.seq, -derep.idx, -asv.seq, -asv.idx)
}

#' Test if all characters in a character vector are members of an alphabet.
#'
#' @param seq (\code{character}) character string(s) to test
#' @param alphabet (\code{character} with all elements of width 1)
#'
#' @return TRUE if all characters in \code{seq} are also in \code{alphabet}
#' @details This function internally uses regular expressions, so
#' \code{alphabet} should not begin with "^" or contain "\\".  "-", which
#' commonly represents a gap, is handled correctly.
#' @export
has_alphabet <- function(seq, alphabet) {
   regex <- paste0("^[", paste0(alphabet, collapse = ""), "]+$")
   # make sure '-' is not interpreted as defining a character range
   regex <- sub(x = regex, pattern = "-", replacement = "\\-")
   return(all(grepl(pattern = regex, x = seq, perl = TRUE)))
}

#' Calculate consensus of a cluster of sequences.
#'
#' This algorithm assumes that the sequences "should be" identical except for
#' amplification and sequencing errors.  Its main purpose is to calculate a
#' consensus sequence for an amplicon that is too long to use in DADA2 directly,
#' but which has been clustered based on sequence variant identity in one
#' subregion.
#'
#' @param seq (\code{character} vector or \code{\link[Biostrings]{XStringSet}})
#'  The sequences to calculate a consensus for.
#' @param names (\code{character}) If \code{seq} is a \code{character} vector,
#'  names for the sequences.
#' @param ncpus (\code{integer}) Number of CPUs to use.
#' @param ... passed to methods
#'
#' @details The sequences are first aligned using
#' \code{\link[DECIPHER]{AlignSeqs}}. Sequences which are "outliers" in the
#' alignment are then removed by
#' \code{\link[odseq]{odseq}}. If the input sequences were clustered based on
#' DADA2 sequence variants of a variable region, and the sequences were
#' appropriately quality filtered prior to running \code{\link[dada2]{dada}},
#' then outliers should mostly be chimeras.
#'
#' After outlier removal, sites with greater than 50\% gaps are removed, and
#' the most frequent base (ignoring gaps) is chosen at all other sites. If no
#' base has greater than 50\% representation at a position, then an IUPAC
#' ambiguous base representing at least 50\% of the reads at that position is
#' chosen.
#'
#' @return an \code{\link[Biostrings]{XStringSet}} representing the consensus
#' sequence.
#' @export

cluster_consensus <- function(seq, ..., ncpus = 1) UseMethod("cluster_consensus")
#' @param DNA2RNA (logical) whether to convert \code{seq} from DNA to RNA, and use (calculated) RNA secondary structure in alignments.
#' @rdname cluster_consensus
#' @export
cluster_consensus.character <- function(seq, names, DNA2RNA = TRUE, ..., ncpus = 1) {
   seq <- rlang::set_names(seq, names)
   seq <- stats::na.omit(seq)
   if (has_alphabet(seq, Biostrings::DNA_ALPHABET)) {
      seq <- Biostrings::DNAStringSet(seq)
      if (DNA2RNA) {
         seq <- Biostrings::RNAStringSet(seq)
      }
   } else if (has_alphabet(seq, Biostrings::RNA_ALPHABET)) {
      seq <- Biostrings::RNAStringSet(seq)
   }
   cluster_consensus.XStringSet(seq, ncpus)
}

#' @rdname cluster_consensus
#' @export
cluster_consensus.XStringSet <- function(seq, ..., ncpus = 1) {

   if (length(seq) < 3) return(NA_character_)

   if (methods::is(seq, "RNAStringSet")) {
      mult_align_class <- Biostrings::RNAMultipleAlignment
      seqset_class <- "RNAStringSet"
   } else if (methods::is(seq, "DNAStringSet")) {
      mult_align_class <- Biostrings::DNAMultipleAlignment
      seqset_class <- "DNAStringSet"
   } else if (methods::is(seq, "AAStringSet")) {
      mult_align_class <- Biostrings::AAMultipleAlignment
      seqset_class <- "AAStringSet"
   } else {
      stop("Unknown sequence class")
   }

   cat("Calculating consensus of", length(seq), "sequences...\n")
   tictoc::tic("total")
   on.exit(tictoc::toc())

   cat(" Aligning...\n")
   tictoc::tic("  alignment")
   aln <- DECIPHER::AlignSeqs(seq, processors = ncpus, verbose = FALSE)
   tictoc::toc()

   cat(" Removing outliers...\n")
   tictoc::tic("  outliers")
   outliers <- odseq::odseq(Biostrings::RNAMultipleAlignment(aln))
   cat("  -removed", sum(outliers), "/", length(outliers),
       "sequences as outliers.\n")
   aln <- aln[!outliers]
   tictoc::toc()

   cat(" Masking gaps...\n")
   tictoc::tic("  masking")
   aln <- aln %>%
      Biostrings::RNAMultipleAlignment() %>%
      Biostrings::maskGaps(min.fraction = 0.5, min.block.width = 1) %>%
      methods::as("RNAStringSet")
   tictoc::toc()

   cat(" Calculating consensus...\n")
   tictoc::tic("  consensus")
   on.exit(tictoc::toc(), add = TRUE)
   DECIPHER::ConsensusSequence(aln,
                               threshold = 0.5,
                               ambiguity = TRUE,
                               ignoreNonBases = TRUE,
                               includeTerminalGaps = FALSE)
}

#' Extract regions from a set of sequences (maybe with qualities)
#'
#' @param seq (\code{character} (a file name) or a \code{\link[ShortRead:ShortRead-class]{ShortRead}} object) the sequences to extract regions from.
#' @param positions (\code{data.frame}) as returned by \code{\link[rITSx]{itsx}} with \code{positions = TRUE} and \code{read_function} set; should have columns \code{$seq} with sequence IDs (matching those in \code{seq}), \code{$region} giving the name of each region, and \code{$start} and \code{$end} giving the start and stop location, if found, of each region.
#' @param region (\code{character}) The region to extract. Should match a value given in \code{positions$region}.
#' @param region2 (\code{character}) If different from \code{region}, then the entire segment beginning at the start of \code{region} and ending at the end of \code{region2} will be extracted.  For instance, to extract the entire ITS region, use \code{region = 'ITS1', region2 = 'ITS2'}.
#' @param outfile (\code{character}) If given, the output will be written to the filename given in fasta or fastq format.  The format is determined by \code{seq}, not by the extension of \code{outfile}.
#' @param ... Passed to methods.
#'
#' @return (\code{object of class \link[ShortRead:ShortRead-class]{ShortRead} or \link[ShortRead:ShortReadQ-class]{ShortReadQ}}) The requested region from each of the input sequences where it was found.
#' @export
#'
extract_region <- function(seq, positions, region, region2 = region, outfile = NULL, ...)
   UseMethod("extract_region")

#' @param qualityType (\code{character}) fastq file quality encoding; see \code{\link[ShortRead]{readFastq}}.
#' @rdname extract_region
#' @export
extract_region.character <- function(seq, positions, region, region2 = region,
                                     outfile = NULL,
                                     qualityType = "FastqQuality", ...) {
   assertthat::assert_that(assertthat::is.string(seq),
                           file.exists(seq))

   if (grepl(seq, pattern = "\\.fastq(\\.gz)?$")) {
      seq <- ShortRead::readFastq(seq, qualityType = qualityType)
   } else if (grepl(seq, pattern = "\\.(fasta|fa|fst)(\\.gz)?$")) {
      seq <- ShortRead::readFasta(seq) %>%
         ShortRead::ShortRead(sread = seq,
                              id = Biostrings::BStringSet(names(seq)))

   }

   extract_region.ShortRead(seq = seq,
                            positions = positions,
                            region = region,
                            region2 = region2,
                            outfile = outfile,
                            ...)
}

#' @rdname extract_region
#' @export
extract_region.ShortRead <- function(seq, positions, region, region2 = region,
                                     outfile = NULL, ...) {

   assertthat::assert_that(
      assertthat::is.string(region),
      assertthat::is.string(region2),
      assertthat::has_name(positions, "seq"),
      assertthat::has_name(positions, "region"),
      assertthat::has_name(positions, "start"),
      assertthat::has_name(positions, "end"),
      is.character(positions$seq),
      is.character(positions$region),
      is.integer(positions$start),
      is.integer(positions$end))

   if (!is.null(outfile)) {
      assertthat::assert_that(assertthat::is.string(outfile))
      #create the output directory if needed
      dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
      assertthat::assert_that(dir.exists(dirname(outfile)))

      # make sure the file exists even if we don't have anything to write.
      if (file.exists(outfile)) file.remove(outfile)
      if (methods::is(seq, "ShortReadQ")) {
         ShortRead::writeFastq(ShortRead::ShortReadQ(), outfile)
      } else {
         ShortRead::writeFasta(ShortRead::ShortRead())
      }
   }

   # if the region is "full", then we don't need to cut anything.
   if (region %in% c("full", "long", "short")) {
      if (!is.null(outfile)) {
         if (methods::is(seq, "ShortReadQ")) {
            ShortRead::writeFastq(seq, outfile, mode = "a")
         } else {
            ShortRead::writeFasta(seq, outfile, mode = "a")
         }
         return(seq)
      }
   }

   p <- positions %>%
      tidyr::gather(
         key = "border", value = "loc",
         start, end) %>%
      dplyr::filter((border == "start" & region == !!region) |
                       (border == "end" & region == region2)) %>%
      dplyr::select(-region) %>%
      tidyr::spread(key = "border", value = "loc") %>%
      dplyr::filter(!is.na(start),
                    start > 0,
                    !is.na(end),
                    end > 0,
                    # end <= readr::parse_number(length),
                    end > start)

   idx <- tibble::tibble(seq = as.character(seq@id),
                 idx = seq_along(seq)) %>%
      dplyr::left_join(dplyr::select(p, "seq"), ., by = "seq") %>%
      dplyr::pull(idx)

   if (nrow(p)) {
      out <- ShortRead::narrow(seq[idx], start = p$start, end = p$end)
      if (!is.null(outfile)) {
         if (methods::is(out, "ShortReadQ")) {
            ShortRead::writeFastq(out, outfile, mode = "a")
         } else {
            ShortRead::writeFasta(out, outfile, mode = "a")
         }
      }
   } else {
      out <- seq[FALSE]
   }
   return(out)
}

# TODO add methods for DNAStringSet, QualityScaledDNAStringSet