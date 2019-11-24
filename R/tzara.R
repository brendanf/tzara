#' @import utils
utils::globalVariables(c("."))
#' @importFrom rlang .data
#' @importFrom futile.logger flog.namespace flog.trace flog.debug flog.info
#' @importFrom stringr str_extract str_replace str_c
#' @importFrom assertthat assert_that
#' @importFrom tibble tibble
#' @importFrom magrittr %>%
`%>%`


.onLoad <- function(libname, pkgname) { #nolint
   backports::import(pkgname)
}


# Combine futile.logger and tictoc by putting the output of toc in a log message
flog_toc <- function(
   level = c("INFO", "TRACE", "DEBUG", "WARN", "ERROR",
             "FATAL", "CARP"),
   func.toc = tictoc::toc.outmsg, #nolint
   ...
) {
   level <- match.arg(level)
   ffunc <-
      switch(level,
             TRACE = futile.logger::flog.trace,
             DEBUG = futile.logger::flog.debug,
             INFO = futile.logger::flog.info,
             WARN = futile.logger::flog.warn,
             ERROR = futile.logger::flog.error,
             FATAL = futile.logger::flog.fatal,
             CARP = futile.logger::flog.carp
      )
   toc <- tictoc::toc(quiet = TRUE)
   ffunc(func.toc(toc$tic, toc$toc, toc$msg), ...)
}

#' Combine DADA2 \code{\link[dada2:derep-class]{derep}} objects into master map
#'
#' @param dereps a (possibly named) \code{list} of
#'        (\code{\link[dada2:derep-class]{derep}} objects), or a
#'        \code{\link[tibble]{tibble}} with a column "derep" containing such a
#'        list.
#' @param .data a \code{\link[tibble]{tibble}} with the same number of rows as
#'        the length of \code{dereps}.
#' @param ... additional columns to add to the output.
#'
#' @details To be useful for further analysis, each sequence should be uniquely
#'  identified.  This can be done in several ways:
#'   \itemize{
#'     \item{\code{dereps} is a named \code{list} with unique names;}
#'     \item{\code{dereps} is a \code{\link[tibble]{tibble}} with columns (other
#'           than "derep") which uniquely identify the rows;}
#'     \item{\code{.data} is provided and its rows are unique; or}
#'     \item{\code{...} is provided and its combinations are unique.}}
#'
#' @return \code{list} with two members: \describe{
#'   \item{\code{$map} (\code{\link[tibble]{tibble}})}{with columns:
#'        "\code{file}" (\code{character}), "\code{idx}" (\code{integer}), and
#'        "\code{map}" (\code{integer}), giving the mapping from the
#'        "\code{idx}"th sequence in "\code{file}" to a sequence in
#'        "\code{fasta}"}
#'   \item{\code{$fasta} (\code{\link[Biostrings]{DNAStringSet}})}{all unique
#'        sequences; the name of each sequence is an \code{integer} which
#'        matches a value in \code{map$newmap}}
#'   }
#' @export
combine_derep <- function(dereps, .data = NULL, ...) {

   # handle different input types to get a tibble with the derep objects in one
   # column
   if (is.data.frame(dereps)) {
      if (length(list(...)) > 0) {
         dereps <- dplyr::bind_cols(
            tibble(.placeholder = seq_along(nrow(dereps)),
                        ...),
            dereps)
         dereps <- dplyr::select(dereps, -".placeholder")
      }
   } else {
      n <- names(dereps)
      dereps <- tibble(..., derep = dereps)
      if (!hasName(dereps, "name") && !is.null(n)) dereps$name <- n
   }

   if (!missing(.data)) {
      dereps <- dplyr::bind_cols(.data, dereps)
   }

   # Check that our derep objects are uniquely identified
   gps <- setdiff(names(dereps), "derep")
   assert_that(dplyr::n_distinct(dereps[gps]) == nrow(dereps))

   # get all the old mappings
   # preserve the sequence names as "seq.id" if they are present
   oldmap <- dereps
   oldmap[["oldmap"]] <- lapply(oldmap[["derep"]], `[[`, "map")
   nestedcols <- "oldmap"
   if (all(vapply(oldmap[["derep"]], assertthat::has_name, TRUE, "names"))) {
      oldmap[["seq.id"]] <- lapply(oldmap[["derep"]], `[[`, "names")
      nestedcols <- c(nestedcols, "seq.id")
   }
   oldmap <- dplyr::select(oldmap, -"derep") %>%
      tidyr::unnest(cols = nestedcols) %>%
      dplyr::group_by_at(gps) %>%
      dplyr::mutate(idx = 1:dplyr::n()) %>%
      dplyr::ungroup()

   # get the old unique sequences
   olduniques <- dereps %>%
      dplyr::mutate_at("derep",
                       ~purrr::map(., .f = ~tibble(seq = names(.$uniques),
                                                           n = .$uniques))) %>%
      tidyr::unnest(cols = "derep")

   # combine duplicate sequences among all files.
   newuniques <- olduniques %>%
      dplyr::group_by(seq) %>%
      dplyr::summarize(n = sum(n)) %>%
      dplyr::arrange(dplyr::desc(n)) %>%
      dplyr::transmute(seq = seq,
                       newmap = seq_along(.data$seq))

   # create a mapping from the unique sequences list in each file
   # to the master unique sequence list
   newderep <- olduniques %>% {
      dplyr::left_join(
      dplyr::select(., dplyr::one_of(gps), seq) %>%
         dplyr::group_by_at(gps) %>%
         dplyr::mutate(oldmap = seq_along(seq)) %>%
         dplyr::ungroup(),
      newuniques,
      by = "seq")
   }

   out <- list()
   # map from the individual sequences in each file to the master unique
   # sequence list
   out$map <- dplyr::left_join(oldmap,
                               dplyr::select(newderep, dplyr::one_of(gps),
                                             "oldmap", "newmap"),
                               by = c(gps, "oldmap")) %>%
      dplyr::select(-oldmap, map = "newmap")
   #unique sequence list
   out$fasta <- Biostrings::DNAStringSet(
      x = purrr::set_names(newuniques$seq, newuniques$newmap)
   )
   class(out) <- c("multiderep", class(out))
   return(out)
}


#' Produce a map between denoised sequences and their original identifiers.
#'
#' @param derep a \code{\link[dada2]{derep-class}} object or list of such
#'        objects
#' @param dada (\link[dada2]{dada-class} object or list of such objects) the
#'        results of a call to \code{\link[dada2]{dada}} on \code{derep}
#' @param ... additional columns to add to the output.  Names included in the
#'        output by default should be avoided.
#'
#' @details Columns \code{$derep.seq} and \code{$dada.seq} contain one sequence
#' per read as plain character strings. Keeping a separate, dereplicated list of
#' sequences and storing references to them may seem like it would be more
#' memory efficient, but it is not necessary to do this explicitly, because this
#' is what R already does with \code{character} vectors; only one copy of each
#' unique string is actually kept in memory, and everything else is pointers.
#'
#' @return a \code{\link[tibble]{tibble}} with columns:
#'   \describe{
#'     \item{\code{$name} (\code{character})}{names from \code{derep}, if it is
#'     a named list of \code{\link[=add_derep_names]{named_derep}} objects.
#'     Otherwise absent unless provided in \code{...}.}
#'     \item{\code{...} (\code{character})}{any additional arguments, passed on
#'       to \code{\link[tibble]{tibble}}}
#'     \item{\code{$seq.id} (\code{character})}{the sequence identifiers from
#'       the original fasta/q file.}
#'     \item{\code{$derep.idx} (\code{integer})}{the index of each sequence in
#'       \code{derep$uniques}.}
#'     \item{\code{$derep.seq} (\code{character})}{the sequences.}
#'     \item{\code{$dada.seq} (\code{character})}{the denoised sequences.}
#' }
#' @export

dadamap <- function(derep, dada, ...) UseMethod("dadamap")
#' @export
dadamap.derep <- function(derep, dada, ...) {
   m <- tibble(
      seq.id = derep$names,
      derep.idx = derep$map,
      derep.seq = names(derep$uniques)[.data$derep.idx],
      ...
   )

   m <- dplyr::left_join(
      m,
      tibble(
         dada.idx = dada$map,
         derep.idx = seq_along(.data$dada.idx),
         dada.seq = dada$sequence[.data$dada.idx]
      ),
      by = "derep.idx"
   )
   class(m) <- c("dadamap", class(m))
   m
}

#' @export
dadamap.list <- function(derep, dada, ...) {
   assert_that(assertthat::are_equal(length(derep), length(dada)))
   assert_that(all(
      purrr::map_lgl(derep, is.null) == purrr::map_lgl(dada, is.null)
   ))
   for (i in seq_along(derep)) {
      assert_that(methods::is(derep[[i]], "derep") || is.null(derep[[i]]))
      assert_that(methods::is(dada[[i]], "dada") || is.null(dada[[i]]))
   }

   args <- list(..., derep = derep, dada = dada)
   if (!is.null(names(derep))) {
      args[["name"]] <- names(derep)
      args <- args[tidyselect::vars_select(names(args), "name",
                                           tidyselect::everything())]
   }
   args <- do.call(tibble::tibble, args) %>%
      dplyr::filter(!purrr::map_lgl(derep, is.null))

   out <- purrr::pmap_dfr(args, dadamap.derep)
   class(out) <- c("dadamap", class(out))
   out
}


#' Individually hash biological sequences
#'
#' @param seq (\code{character} or \code{\link[Biostrings]{XStringSet}}) the
#'        sequences to hash.
#' @param algo (\code{character}) a hash algorithm supported by
#'        \code{\link[digest]{digest}}. default: "xxhash32"
#' @param len (\code{integer}) number of characters to keep from each hash
#'        string. \code{NA} (the default) to keep all characters.
#' @param preserve_na (\code{logical}) If \code{TRUE}, \code{NA} values in
#'        \code{seq} are preserved as \code{NA} in the output.  If
#'        \code{FALSE}, then \code{NA} is passed to
#'        \code{\link[digest]{digest}}, which results in a valid hash.
#'
#' @return a \code{character} vector of the same length as \code{seq},
#'         with the hashed sequences.
#' @export
seqhash <- function(seq, algo = "xxhash32", len = NA, preserve_na = TRUE) {
   UseMethod("seqhash")
}

#' @rdname seqhash
#' @export
seqhash.character <- function(
   seq,
   algo = "xxhash32",
   len = NA,
   preserve_na = TRUE
) {
   h <- vapply(seq, digest::digest, "", algo = algo)
   if (preserve_na) h[is.na(seq)] <- NA_character_
   if (is.na(len)) {
      return(h)
   } else {
      substring(h, 1, len)
   }
}

#' @rdname seqhash
#' @export
seqhash.XStringSet <- function(
   seq,
   algo = "xxhash32",
   len = NA,
   preserve_na = TRUE
) {
   seqhash.character(as.character(seq), algo = algo, len = len, preserve_na)
}

#' Add sequence names to a derep object
#'
#' @param derep (object of class \code{\link[dada2:derep-class]{derep}} or a
#'        \code{list} of  such objects) object(s) to add names to.
#' @param ... passed to methods
#'
#' @return (object of class \code{derep}, or a list of such objects) a
#'         shallow copy of \code{derep}, with an additional member
#'         "\code{$names}", giving the identifiers for the sequences from the
#'         original fasta/q file.
#' @export
add_derep_names <- function(derep, ...) UseMethod("add_derep_names")

#' @rdname add_derep_names
#' @param filename (\code{character}) name of fasta/q file that the
#'     \code{\link[dada2:derep-class]{derep}} object was derived from.
#' @export
add_derep_names.derep <- function(derep, filename, ...) {
   assert_that(file.exists(filename))
   fqs <- ShortRead::FastqStreamer(filename, n = 1e4)
   on.exit(close(fqs))
   seqids <- character(0)
   while (length(fq <- ShortRead::yield(fqs, qualityType = "FastqQuality"))) {
      seqids <- c(seqids, as.character(fq@id))
   }
   derep[["names"]] <- seqids
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
   assert_that(assertthat::are_equal(length(derep), length(names)))
   derep2 <- mapply(add_derep_names.derep, derep, filenames)
   names(derep2) <- filenames
   return(derep2)
}

#' Summarize one or more reads as a \code{\link[tibble]{tibble}}.
#'
#' @param sread (\code{\link[ShortRead:ShortReadQ-class]{ShortReadQ}} object or
#'        list of such objects) as produced by
#'        \code{\link[ShortRead]{readFastq}}.
#' @param max_ee (\code{numeric}) filter out reads with expected error greater
#'        than \code{max_ee}.  Default: Inf (no filtering)
#' @param ... (any vector) additional columns to add to the output
#'        \code{\link[tibble]{tibble}}.  These should have length 1 or, if
#'        \code{sread} is a named list of
#'        \code{\link[ShortRead:ShortReadQ-class]{ShortReadQ}}, the length of
#'        \code{sread}.  Recycling behavior is as in
#'        \code{\link[tibble]{tibble}}. Avoid the name "\code{name}" if
#'        \code{sread} is a named \code{list}
#'
#' @return a \code{\link[tibble]{tibble}} with columns: \describe{
#'    \item{\code{$seq.id} (\code{character})}{sequence IDs, typically from the
#'         input fastq}
#'    \item{\code{$seq} (\code{character})}{the actual sequence reads}
#'    \item{\code{$name} (\code{character})}{the name of the source
#'         \code{\link[ShortRead:ShortReadQ-class]{ShortReadQ}}, if \code{sread}
#'         was a list. Otherwise absent.}
#'    \item{\code{...}}{any other arguments passed to summarize_sread.}}
#' If the input was empty, not a valid
#' \code{\link[ShortRead:ShortReadQ-class]{ShortReadQ}}, or no sequences passed
#' filtering, a \code{\link[tibble]{tibble}} with zero rows is returned.
#' @export
summarize_sread <- function(sread, ..., max_ee) UseMethod("summarize_sread")

#' @export
summarize_sread.ShortReadQ <- function(sread, ..., max_ee = Inf) {
   if (!methods::is(sread, "ShortReadQ")) {
      return(
         tibble(
            seq.id = character(),
            seq = character()
         )
      )
   }

   out <- tibble(
      seq.id = as.character(sread@id),
      seq = as.character(sread@sread),
      ...)
   ee <- rowSums(10 ^ (-1 * (methods::as(sread@quality, "matrix") / 10)),
                 na.rm = TRUE)
   out <- out[ee <= max_ee, , drop = FALSE]
}

#' @export
summarize_sread.list <- function(sread, ..., max_ee = Inf) {
   out <- tibble(
      .seqs = lapply(sread, summarize_sread.ShortReadQ, max_ee = max_ee),
      ...
      )
   if (!is.null(names(sread)) && !hasName(out, "name")) {
      out$name <- names(sread)
   }
   tidyr::unnest(out, ".seqs")
}

#' Test if all characters in a character vector are members of an alphabet.
#'
#' @param seq (\code{character}) character string(s) to test
#' @param alphabet (\code{character} with all elements of width 1)
#'
#' @return TRUE if all characters in \code{seq} are also in \code{alphabet}
#' @details This function internally uses regular expressions, so
#'          \code{alphabet} should not begin with "^" or contain "\\".  "-",
#'          which commonly represents a gap, is handled correctly.
#' @export
has_alphabet <- function(seq, alphabet) {
   regex <- paste0("^[", paste0(alphabet, collapse = ""), "]+$")
   # make sure '-' is not interpreted as defining a character range
   regex <- sub(x = regex, pattern = "-", replacement = "\\\\-")
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
#'        The sequences to calculate a consensus for.
#' @param nread (\code{integer} vector) For the purposes of calculating the
#'        consensus, consider each read to occur \code{nread} times.  Supplying
#'        unique values for \code{seq} along with the corresponding \code{nread}
#'        is much faster than supplying duplicate reads to
#'        \code{cluster_consensus}.
#' @param names (\code{character}) If \code{seq} is a \code{character} vector,
#'        names for the sequences.
#' @param ncpus (\code{integer}) Number of CPUs to use.
#' @param simplify (\code{logical}) If \code{TRUE}, return an object of the same
#'        type as \code{seq} containing a single sequence representing the
#'        consensus. If \code{FALSE}, an object of the same type as \code{seq}
#'        representing the consensus sequence for reads which were included in
#'        the consensus, or \code{NA_character_} for reads which were initially
#'        \code{NA} or which were removed from the consensus alignment as
#'        outliers.  For the \code{\link[Biostrings]{XStringSet}} method, which
#'        does not allow \code{NA} entries, these elements are missing from the
#'        set (this can be deduced by the names).
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
#' the most frequent letter (ignoring gaps) is chosen at all other sites. If no
#' letter has greater than 50\% representation at a position, then an IUPAC
#' ambiguous base representing at least 50\% of the reads at that position is
#' chosen for nucleotide sequences, or \code{"X"} for amino acids.
#'
#' @return an \code{\link[Biostrings]{XStringSet}} representing the consensus
#' sequence.
#'
#' @export

cluster_consensus <- function(seq, nread = 1, ..., ncpus = 1, simplify = TRUE) {
   UseMethod("cluster_consensus")
}
#' @param dna2rna (logical) whether to convert \code{seq} from DNA to RNA, and
#'        use (calculated) RNA secondary structure in alignments.
#' @rdname cluster_consensus
#' @export
cluster_consensus.character <- function(seq, nread = 1, names = names(seq), dna2rna = TRUE,
                                        ..., ncpus = 1, simplify = TRUE) {
   seq <- rlang::set_names(seq, names)
   xss <- seq[!is.na(seq)]
   nread <- nread[!is.na(seq)]
   if (has_alphabet(xss, Biostrings::DNA_ALPHABET)) {
      xss <- Biostrings::DNAStringSet(xss)
      if (dna2rna) {
         xss <- Biostrings::RNAStringSet(xss)
      }
   } else if (has_alphabet(xss, Biostrings::RNA_ALPHABET)) {
      xss <- Biostrings::RNAStringSet(xss)
   }
   result <-
      cluster_consensus.XStringSet(xss, nread = nread, ncpus = ncpus, simplify = simplify)
   if (simplify) {
      if (length(result) == 1) return(as.character(result))
      return(NA_character_)
   }
   seq[] <- NA_character_
   seq[names(result)] <- as.character(result)
   seq
}

#' @rdname cluster_consensus
#' @export
cluster_consensus.XStringSet <- function(seq, nread = 1, ..., ncpus = 1, simplify = TRUE) {

   assertthat::assert_that(length(nread) == 1 | length(nread) == length(seq))
   if (sum(nread) < 3) return(seq[FALSE])
   if (length(seq) == 1) return(seq)

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

   flog.info("Calculating consensus of %d sequences...", length(seq))
   tictoc::tic("cluster_consensus")
   on.exit(flog_toc("DEBUG"))

   flog.debug("Aligning...")
   tictoc::tic("cluster_consensus:aligning")
   aln <- DECIPHER::AlignSeqs(seq, processors = ncpus, verbose = FALSE)
   flog_toc("TRACE")

   flog.debug("Removing outliers...")
   tictoc::tic("cluster_consensus:outliers")
   outliers <- odseq(mult_align_class(aln), weights = nread)
   flog.trace("Removed %d/%d sequences as outliers.",
                             sum(outliers), length(outliers))
   aln <- aln[!outliers]
   nread <- nread[!outliers]
   flog_toc("TRACE")

   flog.trace("Masking gaps...")
   tictoc::tic("cluster_consensus:masking")
   aln2 <- rep(aln, nread)
   aln2 <- aln2 %>%
      mult_align_class() %>%
      Biostrings::maskGaps(min.fraction = 0.5, min.block.width = 1) %>%
      methods::as(seqset_class)
   flog_toc("TRACE")

   flog.trace("Calculating consensus...")
   tictoc::tic("cluster_consensus:consensus")
   on.exit(flog_toc("TRACE"), add = TRUE, after = FALSE)
   result <- DECIPHER::ConsensusSequence(
      aln2,
      threshold = 0.5,
      ambiguity = TRUE,
      ignoreNonBases = TRUE,
      includeTerminalGaps = FALSE
   )
   if (simplify) return(result)

   result <- rep(result, length(aln))
   names(result) <- names(aln)
   result
}

#' Extract regions from a set of sequences (maybe with qualities)
#'
#' @param seq (\code{character} (a file name) or a
#'        \code{\link[ShortRead:ShortRead-class]{ShortRead}} object) the
#'        sequences to extract regions from.
#' @param positions (\code{data.frame}) as returned by \code{\link[rITSx]{itsx}}
#'        with \code{positions = TRUE} and \code{read_function} set; should have
#'        columns \code{$seq} with sequence IDs (matching those in \code{seq}),
#'        \code{$region} giving the name of each region, and \code{$start} and
#'        \code{$end} giving the start and stop location, if found, of each
#'        region.
#' @param region (\code{character}) The region to extract. Should match a value
#'        given in \code{positions$region}.
#' @param region2 (\code{character}) If different from \code{region}, then the
#'        entire segment beginning at the start of \code{region} and ending at
#'        the end of \code{region2} will be extracted.  For instance, to extract
#'        the entire ITS region, use \code{region = 'ITS1', region2 = 'ITS2'}.
#' @param outfile (\code{character}) If given, the output will be written to the
#'        filename given in fasta or fastq format.  The format is determined by
#'        \code{seq}, not by the extension of \code{outfile}.
#' @param ... Passed to methods.
#'
#' @return (\code{object of class \link[ShortRead:ShortRead-class]{ShortRead} or
#'         \link[ShortRead:ShortReadQ-class]{ShortReadQ}}) The requested region
#'         from each of the input sequences where it was found.
#' @export
#'
extract_region <- function(seq, positions, region, region2 = region,
                           outfile = NULL, ...)
   UseMethod("extract_region")

#' @param qualityType (\code{character} scalar) fastq file quality encoding; see
#'        \code{\link[ShortRead]{readFastq}}.
#' @param append (\code{logical} scalar) if \code{TRUE}, then data is appended
#'        to \code{outfile}; if \code{FALSE}, existing data in \code{outfile} is
#'        overwritten.
#' @rdname extract_region
#' @export
extract_region.character <- function(seq, positions, region, region2 = region,
                                     outfile = NULL,
                                     qualityType = "FastqQuality", #nolint
                                     append = FALSE, ...) {
   assert_that(assertthat::is.flag(append))

   if (length(seq) > 1) {
      assert_that(length(positions) == length(seq))
      if (!is.null(outfile) && !append) unlink(outfile)
      out <- purrr::map2(
         .x = seq,
         .y = positions,
         .f = extract_region.character,
         region = region,
         region2 = region2,
         outfile = outfile,
         qualityType = qualityType,
         append = TRUE,
         ...
      )
      return(purrr::reduce(out, ShortRead::append))
   }

   assert_that(assertthat::is.string(seq),
                           file.exists(seq))

   assert_that(is.data.frame(positions) || is.list(positions))
   if (!is.data.frame(positions)) {
      assert_that(length(positions) == length(seq))
      positions <- positions[[1]]
      assert_that(is.data.frame(positions))
   }

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
                            append = append,
                            ...)
}

#' @rdname extract_region
#' @export
extract_region.ShortRead <- function(seq, positions, region, region2 = region,
                                     outfile = NULL, append = FALSE, ...) {

   assert_that(
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
      assert_that(assertthat::is.string(outfile))
      #create the output directory if needed
      dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
      assert_that(dir.exists(dirname(outfile)))

      # make sure the file exists even if we don't have anything to write.
      if (file.exists(outfile) && !isTRUE(append)) file.remove(outfile)
      if (!file.exists(outfile)) {
         if (methods::is(seq, "ShortReadQ")) {
            ShortRead::writeFastq(ShortRead::ShortReadQ(), outfile)
         } else {
            ShortRead::writeFasta(ShortRead::ShortRead(), outfile)
         }
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
      }
      return(seq)
   }

   p <- positions %>%
      tidyr::gather(
         key = "border", value = "loc",
         "start", "end") %>%
      dplyr::filter(
         (.data$border == "start" & .data$region == !!region) |
            (.data$border == "end" & .data$region == region2)
      ) %>%
      dplyr::select(-"region") %>%
      tidyr::spread(key = "border", value = "loc") %>%
      dplyr::filter(
         !is.na(.data$start),
         .data$start > 0,
         !is.na(.data$end),
         .data$end > 0,
         # end <= readr::parse_number(length),
         .data$end > .data$start
      )

   idx <- tibble(
      seq = as.character(seq@id),
      idx = seq_along(.data$seq)
   ) %>%
      dplyr::left_join(dplyr::select(p, "seq"), ., by = "seq") %>%
      dplyr::pull("idx")

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

#' @rdname extract_region
#' @export
extract_region.list <- function(seq, positions, region, region2 = region,
                                     outfile = NULL, ...) {

   assert_that(length(positions) == length(seq))
   if (!is.null(outfile) && !append) unlink(outfile)
   out <- purrr::map2(
      .x = seq,
      .y = positions,
      .f = extract_region,
      region = region,
      region2 = region2,
      outfile = outfile,
      append = TRUE,
      ...
   )
   purrr::reduce(out, ShortRead::append)
}

str_modify <- function(x, regex, replace, ...) {
   if (!is.null(regex) && !is.na(regex)) {
      assert_that(assertthat::is.string(regex))
      if (is.null(replace) || is.na(replace)) {
         str_extract(x, regex, ...)
      } else {
         assert_that(assertthat::is.string(replace))
         str_replace(x, regex, replace, ...)
      }
   }
}

# TODO add methods for DNAStringSet, QualityScaledDNAStringSet

#' Reconstruct a longer region out of ASVs or consensus sequence of individual
#' domains.
#'
#' The sequences from each denoised sub-region/domain are concatenated to create
#' a denoised sequence
#' for the long region.  Additionally, de-novo bimera detection is performed
#' using \code{\link[dada2]{isBimeraDenovo}} or
#' \code{\link[dada2]{isBimeraDenovoTable}} on
#' sets of three consecutive sub-regions/domains; in the intended application,
#' these sets will be variable--conserved--variable.
#'
#' When not all sub-regions/domains for a given read have been successfully
#' denoised with DADA, then the missing regions are constructed using
#' \code{\link{cluster_consensus}}.
#'
#' @param seqtabs (\code{list} of \code{data.frame}) with columns
#'  \code{read_column}, \code{asv_column}, and optionally \code{sample_column}.
#'  Any additional columns are ignored.  \code{read_column} should give a unique
#'  ID for each sequencing read, and \code{asv_column} should give the denoised
#'  sequence for the read.
#' @param regions  (\code{character} vector with the same length as
#'  \code{seqtabs}) The names of the regions/domains represented by each of the
#'  tables in \code{seqtabs}.  If not supplied, then \code{seqtabs} should be
#'  named by the regions.
#' @param regions_regex (\code{character} scalar, or \code{NULL})
#'  A \link[stringi:stringi-search-regex]{regular expression}. If
#'  \code{regions_regex} is given but \code{regions_replace} is not, then only
#'  the part of the entries in \code{regions} matching the regex
#'  are used to define samples (using \code{\link[stringr]{str_extract}}).  If
#'  \code{regions_replace} is also used, then the regex is instead replaced by
#'  \code{regions_replace} (using \code{\link[stringr]{str_replace}}).
#'  \code{NA_character} is treated the same way as \code{NULL}.
#' @param regions_replace (\code{character} scalar, or \code{NULL})
#'  Replacement string for \code{regions_regex}.
#'  \code{NA_character} is treated the same way as \code{NULL}.
#' @param output (\code{character} scalar or named list of \code{character}
#'  vectors) If a \code{character} scalar, then the name to be used for the
#'  (single) output region.  In this case the region will be the concatenation
#'  of all the regions in \code{order}.  Alternatively, a list where the names
#'  are the names of the output regions, and the values are \code{character}
#'  vectors giving the regions which should be concatenated for each output
#'  region.
#' @param use_output (one of \code{"first"}, \code{"second"}, or \code{"no"}) If
#'        one of the regions given by \code{output} is also present in
#'        \code{seqtabs}, then the \code{seqtabs} version is used preferentially
#'        \code{use_output == "first"}, as a backup value when one of the
#'        subregions/domains is missing if \code{use_output == "second"}, or not
#'        at all if \code{use_output == "no"}.
#' @param order (\code{character} vector) The order in which the
#'        sub-regions/domains should be concatenated to produce the output(s).
#' @param read_column (\code{character} scalar) Column name from the
#'        \code{seqtabs} which uniquely identifies each read (but different
#'        regions extracted from the same read should have the same ID.)
#' @param asv_column (\code{character} scalar) Column name from the
#'        \code{seqtabs} which gives the denoised sequences.
#' @param rawtabs (\code{list} of \code{data.frame}) Data sources of the same
#'        format as \code{seqtabs}, with columns \code{read_column} and
#'        \code{raw_column}.  These should be of the same number as
#'        \code{seqtabs}, and correspond to the sub-regions/domains specified in
#'        \code{regions}.  The default is to look for \code{raw_column} in
#'        \code{seqtabs}.
#' @param raw_column (\code{character} scalar, or \code{NULL})
#'        Column name from the \code{seqtabs} which gives the raw sequences.  If
#'        \code{NULL} or \code{NA_character_}, then consensus sequences will not
#'        be used as a  backup when no denoised sequence is present.
#' @param raw_regions  (\code{character} vector with the same length as
#'        \code{rawtabs}) The names of the regions/domains represented by each
#'        of the tables in \code{rawtabs}.  These will be processed using
#'        \code{regions_regex} and \code{regions_replace}, if given.
#' @param sample_column (\code{character} scalar, or \code{NULL}) An optional
#'        column name from the \code{seqtabs} which identifies which sample each
#'        sequence is from.  If given, this is used (after possible modification
#'        by \code{sample_regex} and \code{sample_replace}) to identify
#'        different samples for \code{\link[dada2]{isBimeraDenovoTable}}.
#'        \code{NA_character} is treated the same way as \code{NULL}.
#' @param sample_regex (\code{character} scalar, or \code{NULL}) A
#'        \link[stringi:stringi-search-regex]{regular expression}. If
#'        \code{sample_regex} is given but \code{sample_replace} is not, then
#'        only the part of the entries in \code{sample_column} matching the
#'        regex are used to define samples (using
#'        \code{\link[stringr]{str_extract}}).  If \code{sample_replace} is also
#'        used, then the regex is instead replaced by \code{sample_replace}
#'        (using \code{\link[stringr]{str_replace}}). \code{NA_character} is
#'        treated the same way as \code{NULL}.
#' @param sample_replace (\code{character} scalar, or \code{NULL})
#'        Replacement string for \code{sample_regex}. \code{NA_character} is
#'        treated the same way as \code{NULL}.
#' @param chimera_offset (\code{integer}) By default, bimeras are checked for
#'        sub-region/domains 1, 2, 3; 3, 4, 5; 5, 6, 7; etc. This is appropriate
#'        if the domains alternate variable, conserved, variable, etc.  If a
#'        more conserved domain is first, use \code{chimera_offset = 1}.
#' @param allow_map (\code{logical} scalar) If \code{TRUE} and if \code{asvs}
#'        contains non-missing values, attempt to map each raw read without a
#'        corresponding ASV to the nearest ASV.
#' @param allow_consensus (\code{logical} scalar) If \code{TRUE} and if
#'        \code{allow_map} is \code{FALSE} or there are no non-missing values in
#'        \code{asvs}, then attempt to make a consensus of all raw reads.
#' @param allow_raw (\code{logical} scalar) If \code{TRUE}, then after mapping
#'        and/or consensus building, remaining raw reads are taken as they are.
#'        If \code{FALSE}, the corresponding results will be \code{NA}.
#' @param ... additional arguments passed to \code{\link[dada2]{isBimeraDenovo}}
#'        or \code{\link[dada2]{isBimeraDenovoTable}}.
#'
#' @return a \code{\link[tibble]{tibble}} with column "\code{seq.id}" and
#'        \code{sample_column} (if given), as well as one column for each value
#'        of \code{regions} and \code{output}, representing the
#'        sub-regions/domains and the concatenated full region.
#' @export
reconstruct <- function(
   seqtabs,
   regions = names(seqtabs),
   regions_regex = NULL,
   regions_replace = NULL,
   output = "concat",
   use_output = c("first", "second", "no"),
   order = setdiff(regions, output),
   read_column = "seq.id",
   asv_column = "dada.seq",
   rawtabs = seqtabs,
   raw_column = NULL,
   raw_regions = names(rawtabs),
   sample_column = NULL,
   sample_regex = NULL,
   sample_replace = NULL,
   chimera_offset = 0,
   allow_map = TRUE,
   allow_consensus = TRUE,
   allow_raw = FALSE,
   ...
) {
   assert_that(
      is_string_or_missing(raw_column),
      is.character(order)
   )

   use_output <- match.arg(use_output)

   if (is.null(raw_column) || is.na(raw_column)) raw_column <- NULL

   regions <- str_modify(regions, regions_regex, regions_replace)

   assert_that(all(order %in% regions))

   if (is.character(output)) {
      assert_that(assertthat::is.string(output))
      output <- magrittr::set_names(list(order), output)
   } else {
      assert_that(
         is.list(output),
         rlang::is_named(output)
      )
   }
   assert_that(
      all(purrr::map_lgl(output, is.character)),
      all(purrr::map_lgl(output, ~all(. %in% order)))
   )

   region_table <- assemble_region_table(
      seqtabs = seqtabs,
      regions = regions,
      output = output,
      order = order,
      read_column = read_column,
      seq_column = asv_column,
      sample_column = sample_column,
      sample_regex = sample_regex,
      sample_replace = sample_replace
   )

   chims <- find_all_region_chimeras(
      region_table = region_table,
      order = order,
      sample_column = sample_column,
      read_column = read_column,
      chimera_offset = chimera_offset)

   if (!is.null(raw_column)) {
      raw_regions <- str_modify(regions, regions_regex, regions_replace)
      raw_table <- assemble_region_table(
         seqtabs = rawtabs,
         regions = raw_regions,
         output = output,
         order = order,
         read_column = read_column,
         seq_column = raw_column,
         sample_column = sample_column,
         sample_regex = sample_regex,
         sample_replace = sample_replace
      )
      raw_table <- dplyr::semi_join(raw_table, region_table, by = read_column)

      region_table <- consensus_missing_regions(
         region_table = region_table,
         raw_table = raw_table,
         order = order,
         read_column = read_column,
         allow_map = allow_map,
         allow_consensus = allow_consensus,
         allow_raw = allow_raw,
         ...
      )
   }

   for (o in names(output)) {
      this_chims <- purrr::map_lgl(chims$chimset, ~all(. %in% output[[o]]))
      this_chims <- chims$chims[this_chims]
      this_chims <- unlist(this_chims)
      this_chims <- unique(this_chims)
      this_chims <- which(region_table[[read_column]] %in% this_chims)
      region_table <- reconstruct_region(
         region_table = region_table,
         output = o,
         use_output = use_output,
         order = output[[o]],
         chims = this_chims
      )
   }

   region_table
}

reconstruct_region <- function(region_table, output, order, use_output,
                               chims = integer(0)) {
   if (output %in% names(region_table) && use_output != "no") {
      region_table[["_concat_"]] <-
         do.call(
            str_c,
            region_table[, order]
         )
      if (use_output == "first") {
         region_table[[output]] <-
            dplyr::coalesce(
               region_table[[output]],
               region_table[["_concat_"]]
            )
      } else {
         region_table[[output]] <-
            dplyr::coalesce(
               region_table[["_concat_"]],
               region_table[[output]]
            )
      }
      region_table[["_concat_"]] <- NULL
   } else {
      region_table[[output]] <-
         do.call(
            str_c,
            region_table[, order]
         )
   }
   region_table[[output]][chims] <- NA_character_
   region_table
}

#' Combine raw and denoised reads from multiple regions of the same sequences.
#'
#' @param dadamap (\code{\link{dadamap}} object)
#' @param rawdata (\code{\link[tibble]{tibble}}) as returned by
#'        \code{\link{summarize_sread}}
#'
#' @details Both of the inputs should be annotated with identifying columns to
#' uniquely identify the source \code{\link[dada2:derep-class]{derep}} objects;
#' this should happen automatically if \code{\link[dada2]{derepFastq}} and
#' \code{\link[ShortRead]{readFastq}} are called on a list of filenames (there
#' will be a column "name" in the outputs of \code{\link{dadamap}} and
#' \code{\link{summarize_sread}}).
#'
#' @return a \code{\link[tibble]{tibble}} giving the sequences for each region
#'         for each read (uniquely identified by all other columns)
#' @export
combine_bigmaps <- function(dadamap, rawdata) {
   joincols <- intersect(colnames(dadamap), colnames(rawdata))
   dplyr::full_join(dadamap, rawdata, by = joincols) %>%
      dplyr::group_by("seq.id") %>%
      dplyr::filter(any(!is.na(.data$dada.idx))) %>%
      dplyr::mutate(
         seq = dplyr::coalesce(.data$dada.seq, .data$derep.seq, .data$seq)
      ) %>%
      dplyr::select(-"derep.seq", -"derep.idx", -"dada.seq", -"dada.idx")
}

is_string_or_missing <- function(x) {
   is.null(x) || is.na(x) || assertthat::is.string(x)
}

assemble_region_table <- function(
   seqtabs,
   regions = names(seqtabs),
   output,
   order,
   read_column = "seq.id",
   seq_column = "dada.seq",
   sample_column = NULL,
   sample_regex = NULL,
   sample_replace = NULL
) {
   assert_that(
      is.list(output),
      rlang::is_named(output),
      assertthat::is.string(read_column),
      assertthat::is.string(seq_column),
      is_string_or_missing(sample_column),
      is_string_or_missing(sample_regex),
      is_string_or_missing(sample_replace)
   )

   if (is.na(sample_column)) sample_column <- NULL

   seqtabs <- purrr::map2(
      seqtabs,
      regions,
      function(st, reg) magrittr::set_names(
         st[, c(sample_column, read_column, seq_column)],
         c(sample_column, read_column, reg)))
   seqtabs <- purrr::map(
      unique(regions) %>% magrittr::set_names(., .),
      function(r) {
         dplyr::bind_rows(seqtabs[regions == r])
      }
   )

   if (!is.null(sample_column)) {
      seqtabs <- purrr::map(
         seqtabs,
         dplyr::mutate_at,
         sample_column,
         str_modify,
         sample_regex,
         sample_replace
      )
   }

   out <- purrr::reduce(
      seqtabs[order],
      dplyr::full_join,
      by = c(sample_column, read_column)
   )
   for (o in names(output)) {
      if (o %in% regions) {
         out <- dplyr::full_join(
            out,
            seqtabs[[o]],
            by = c(sample_column, read_column)
         )
      }
   }
   out
}

#' Find chimeras in adjacent variable-conserved-variable domain triplets
#'
#' @param region_table (\code{data.frame} containing all the regions in
#'        \code{order})
#' @param order (\code{character} with the regions to check)
#' @param sample_column (\code{character} scalar) extra column in region_table
#'        which identifies different samples; if present, use
#'        \code{\link[dada2]{isBimeraDenovoTable}} instead of
#'        \code{\link[dada2]{isBimeraDenovo}}.
#' @param read_column (\code{character} scalar) column in region_table
#'        which identifies different reads.
#' @param chimera_offset (\code{integer} scalar) use \code{1} if the first
#'        region is conserved.
#' @param ... passed to \code{find_region_chimeras}.
#'
#' @return a \code{\link[tibble]{tibble}} with columns:
#'         \describe{
#'           \item{chimset}{(\code{list} of \code{character}) the set of regions
#'               (usually 3) which were checked in each row.}
#'           \item{chims}{\code{list} of \code{character}) the sequence IDs of
#'               the reads which were identified as chimeric for each set of
#'               regions}}
find_all_region_chimeras <- function(
   region_table,
   order,
   sample_column = NULL,
   read_column = "read_id",
   chimera_offset = 0,
   ...
) {
   result <- tibble(
      first = seq(1 + 0, max(1, length(order) - 2), 2),
      last = pmin(length(order), .data$first + 2),
      chimset = purrr::map2(.data$first, .data$last, ~order[.x:.y]),
      chims = purrr::map(
         .data$chimset,
         find_region_chimeras,
         region_table = region_table,
         sample_column = sample_column,
         read_column = read_column,
         ...
      )
   )
   dplyr::select(result, "chimset", "chims")
}


#' Check for bimeras in a subset of regions
#'
#' @param region_table (\code{data.frame} containing all the regions in
#'        \code{chimset})
#' @param chimset (\code{character} vector, usually of length 3) regions to
#'        concatenate and then check for bimeras.
#' @param sample_column (\code{character} scalar) extra column in region_table
#'        which identifies different samples; if present, use
#'        \code{\link[dada2]{isBimeraDenovoTable}} instead of
#'        \code{\link[dada2]{isBimeraDenovo}}.
#' @param read_column (\code{character} scalar) column in region_table
#'        which identifies different reads.
#' @param ... passed on to \code{\link[dada2]{isBimeraDenovoTable}} or
#'        \code{\link[dada2]{isBimeraDenovo}}.
#'
#' @return an \code{character} vector giving the read IDs of
#'        \code{region_table} which were detected as bimeras.
find_region_chimeras <- function(region_table, chimset, sample_column,
                                 read_column, ...) {
   seqs <- do.call(str_c, region_table[, chimset])
   chimset_name <- paste(chimset, collapse = "--")
   flog.trace("Searching for bimeras in regions %s.", chimset_name)
   tictoc::tic(paste(chimset_name), "bimeras")
   chims <-
      if (!is.null(sample_column)) {
         tibble(sample = region_table[[sample_column]], seq = seqs) %>%
            dplyr::filter(!is.na(seq)) %>%
            dplyr::group_by(sample, seq) %>%
            dplyr::summarize(nread = dplyr::n()) %>%
            dplyr::ungroup() %>%
            tidyr::spread(seq, "nread", fill = 0L) %>%
            tibble::column_to_rownames("sample") %>%
            as.matrix() %>%
            dada2::isBimeraDenovoTable(...)
      } else {
         table(seqs) %>% {
            tibble(abundance = ., sequence = names(.))
         } %>%
            dada2::isBimeraDenovo(...)
      }
   flog_toc("TRACE")
   flog.debug("Found %d/%d bimeric ASVs (%d/%d reads) in regions %s.",
              sum(chims), length(chims),
              sum(seqs %in% names(chims)[chims], na.rm = TRUE),
              sum(!is.na(seqs)),
              chimset_name)
   region_table[[read_column]][seqs %in% names(chims)[chims]]
}

block_consensus <- function(.x, .y, reg, reg2, reg2_raw, read_column, ...) {
   if (length(.y[[reg]]) == 0 || is.na(.y[[reg]])) return(.x)
   .x[[reg2]] <- map_or_consensus(
      asvs = .x[[reg2]],
      raw = .x[[reg2_raw]],
      names = .x[[read_column]],
      ...
   )
   .x
}

consensus_missing_regions <- function(
   region_table,
   raw_table,
   order,
   maxdist = 10,
   read_column = "seq.id",
   ...
) {
   # all combinations of ASV sequences
   combos <- dplyr::group_by_at(region_table, order) %>%
      dplyr::summarize(nreads = dplyr::n()) %>%
      dplyr::ungroup()

   # how many variants do we have for each region?
   # sort decreasing, because we want to start with the most variable region
   region_counts <- vapply(
      order,
      function(x) dplyr::n_distinct(combos[[x]], na.rm = TRUE),
      1L
   ) %>%
      sort(decreasing = TRUE)

   out_table <- dplyr::left_join(region_table, raw_table, by = read_column,
                                 suffix = c("", "_raw"))

   for (reg in names(region_counts)) {
      flog.info("Starting to process %s clusters.", reg)
      tictoc::tic(paste(reg, "clusters"))
      out_table <- dplyr::group_by_at(out_table, reg)
      regstats <- dplyr::summarise_at(
         out_table,
         setdiff(order, reg),
         list(
            n_na = ~sum(is.na(.)),
            all_na = ~all(is.na(.)),
            n = length
         )
      )
      for (reg2 in setdiff(order, reg)) {
         all_na <- paste0(reg2, "_all_na")
         n_na <- paste0(reg2, "_n_na")
         flog.info(
            "Interpolating %s values for %s ASVs.", reg2, reg)
         tictoc::tic(paste0(reg, " clusters (", reg2, ")"))

         flog.info(
            "%d missing %s reads found in %d %s clusters with %s ASVs",
            sum(regstats[[n_na]][!regstats[[all_na]] &
                                    !is.na(regstats[[reg]])]),
            reg2,
            sum(
               !regstats[[all_na]] &
                  regstats[[n_na]] > 0 &
                  !is.na(regstats[[reg]])
            ),
            reg,
            reg2
         )
         flog.info(
            "%d missing %s reads found in %d %s clusters with no %s ASVs.",
            sum(regstats[[n_na]][regstats[[all_na]] & !is.na(regstats[[reg]])]),
            reg2,
            sum(regstats[[all_na]] & !is.na(regstats[[reg]])),
            reg,
            reg2
         )
         reg2_raw <- paste0(reg2, "_raw")
         out_table <- dplyr::group_map(
            .tbl = out_table,
           .f = block_consensus,
            reg = reg,
            reg2 = reg2,
            reg2_raw = reg2_raw,
            read_column = read_column,
            ...,
            keep = TRUE
         ) %>%
            dplyr::bind_rows() %>%
            dplyr::group_by_at(reg)
         regstats <- dplyr::summarise_at(
            out_table,
            setdiff(order, reg),
            list(
               n_na = ~sum(is.na(.)),
               all_na = ~all(is.na(.)),
               n = length
            )
         )
         flog_toc()
         flog.info(
            "%d missing %s reads remain in %d %s clusters with %s ASVs",
            sum(regstats[[n_na]][!regstats[[all_na]] &
                                    !is.na(regstats[[reg]])]),
            reg2,
            sum(
               !regstats[[all_na]] &
                  regstats[[n_na]] > 0 &
                  !is.na(regstats[[reg]])
            ),
            reg,
            reg2
         )
         flog.info(
            "%d missing %s reads remain in %d %s clusters with no %s ASVs.",
            sum(regstats[[n_na]][regstats[[all_na]] & !is.na(regstats[[reg]])]),
            reg2,
            sum(regstats[[all_na]] & !is.na(regstats[[reg]])),
            reg,
            reg2
         )
      }
      flog_toc()
   }
   dplyr::select_at(out_table, names(region_table))
}

#' Replace unmapped raw reads with the nearest ASV
#'
#' @param asvs (\code{character} vector) ASV sequences mapped to a set of reads.
#'    Should be \code{\link{NA_character_}} for reads which did not map to an
#'    ASV.
#' @param raw (\code{character} vector) Raw read sequences for the same set of
#'    reads as \code{asvs}.  May be \code{\link{NA_character_}}
#' @param maxdist (\code{numeric} scalar) Maximum Levenshtein distance between
#'    a raw read and an ASV for the read to be mapped to the ASV.
#'
#' @return a \code{character} vector the same length as \code{asvs}, which has
#'    the closest ASV for each read.
#'
#' @details The value for element \code{i} of the result is determined as
#'  follows:
#'   \enumerate{
#'     \item{\code{asvs[i]} is non-\code{NA}: the value from \code{asvs} is
#'      used.}
#'     \item{\code{asvs[i]} and \code{raw[i]} are both \code{NA}:
#'      \code{\link{NA_character_}}}
#'     \item{\code{asvs[i]} is \code{NA}, \code{raw[i]} is non-\code{NA}:
#'      \enumerate{
#'       \item{\code{raw[i]} is less than \code{maxdist} in edit distance from
#'        at least one of the non-\code{NA} sequences in \code{asvs}: the value
#'        from \code{asvs} which has the smallest edit distance from
#'        \code{raw[i]} is used.}
#'       \item{\code{raw[i]} is not less than \code{maxdist} in edit distance
#'        from at least one of the non-\code{NA} sequences in \code{asvs}:
#'        \code{NA_character_}}}}}
#'
#' @export
map_to_best_asv <- function(asvs, raw, maxdist = 10) {
   assert_that(
      is.character(asvs),
      is.character(raw),
      length(asvs) == length(raw),
      assertthat::is.number(maxdist),
      maxdist >= 0
   )
   if (!anyNA(asvs)) return(asvs)

   unique_asvs <- purrr::discard(unique(asvs), is.na)
   required_raw <- raw[is.na(asvs)]
   unique_raw <- purrr::discard(unique(required_raw), is.na)
   if (length(unique_raw) == 0) return(asvs)

   flog.debug(
      paste("Calculating distance matrix between %d raw sequences",
             "and %d ASVs."),
      length(unique_raw),
      length()
   )
   tictoc::tic("Distance matrix")
   d <- adist(unique_raw, unique_asvs)
   flog_toc("DEBUG")
   replace_table <- tibble(
      raw = unique_raw,
      match = unique_asvs[apply(d, MARGIN = 1, which.min)],
      dist = apply(d, MARGIN = 1, min)
   )
   replace_table <- dplyr::filter(
      replace_table,
      .data$dist <= maxdist
   )
   out <- dplyr::left_join(tibble(raw = raw), replace_table, by = "raw")
   dplyr::coalesce(asvs, out$match)
}


#' Assign consensus sequences to unmapped reads
#'
#' This function is intended to be run on reads from a single sub-region/domain,
#' which have been clustered using a linked sub-region/domain; for instance,
#' ITS1 reads clustered based on identity/similarity of the linked ITS2 reads.
#'
#' If some of the target reads have been mapped to ASVs, then
#' \code{map_or_consensus} attempts to map additional raw reads to the same
#' ASVs using a (potentially) more relaxed criteria than
#' \code{\link[dada2]{dada}}. This is implemented in
#' \code{\link{map_to_best_asv}}.
#'
#' If, on the other hand, none of the input sequences have been assigned to an
#' ASV, then the entire group is taken to represent one cluster, and a consensus
#' sequence for the cluster is determined using \code{\link{cluster_consensus}}.
#' This process will remove outliers (generally chimeric in origin) and assign
#' \code{NA_character_} to the associated reads, as well as any reads which are
#' already \code{NA} due to quality filtering, failed region extraction, etc.
#'
#' @param asvs (\code{character} vector) ASV sequences mapped to a set of reads.
#'    Should be \code{\link{NA_character_}} for reads which did not map to an
#'    ASV.
#' @param raw (\code{character} vector) Raw read sequences for the same set of
#'    reads as \code{asvs}.  May be \code{\link{NA_character_}}
#' @param maxdist (\code{numeric} scalar) Maximum Levenshtein distance between
#'    a raw read and an ASV for the read to be mapped to the ASV.
#' @param allow_map (\code{logical} scalar) If \code{TRUE} and if \code{asvs}
#'    contains non-missing values, attempt to map each raw read without a
#'    corresponding ASV to the nearest ASV.
#' @param allow_consensus (\code{logical} scalar) If \code{TRUE} and if
#'    \code{allow_map} is \code{FALSE} or there are no non-missing values in
#'    \code{asvs}, then attempt to make a consensus of all raw reads.
#' @param allow_raw (\code{logical} scalar) If \code{TRUE}, then after mapping
#'    and/or consensus building, remaining raw reads are taken as they are. If
#'    \code{FALSE}, the corresponding results will be \code{NA}.
#' @param ... passed to \code{\link{cluster_consensus.character}}
#'
#' @return a \code{character} vector the same length as \code{asvs}, which has
#'    the closest ASV for each read if any ASVs are non-missing, or the cluster
#'    consensus values for raw reads which were non-missing and not outliers.
#' @export
map_or_consensus <- function(asvs, raw, maxdist = 10, allow_map = TRUE,
                             allow_consensus = TRUE, allow_raw = FALSE, ...) {
   assert_that(
      is.character(asvs),
      is.character(raw),
      length(asvs) == length(raw),
      assertthat::is.number(maxdist),
      maxdist >= 0,
      assertthat::is.flag(allow_map),
      assertthat::is.flag(allow_consensus),
      assertthat::is.flag(allow_raw)
   )

   result <- asvs
   if (!all(is.na(asvs)) && allow_map) {
      result <- map_to_best_asv(asvs, raw, maxdist)
   } else if (allow_consensus) {
      result <- cluster_consensus(raw, simplify = FALSE, ...)
   }

   if (allow_raw) {
      result <- dplyr::coalesce(result, raw)
   }
   result
}
