# copy of coerce code from LSUx, since it is not required for tzara

#' @importClassesFrom Biostrings DNAStringSet
#' @importClassesFrom ShortRead ShortRead
if (!methods::hasMethod("coerce", c("DNAStringSet", "ShortRead"))) {
   methods::setAs(
      "DNAStringSet",
      "ShortRead",
      function(from) {
         ShortRead::ShortRead(
            magrittr::set_names(from, NULL),
            Biostrings::BStringSet(names(from))
         )
      }
   )
}

if (!methods::hasMethod("coerce", c("ShortRead", "DNAStringSet"))) {
   methods::setAs(
      "ShortRead",
      "DNAStringSet",
      function(from) {
         out <- ShortRead::sread(from)
         names(out) <- as.character(ShortRead::id(from))
         out
      }
   )
}

if (!methods::hasMethod("coerce", c("character", "ShortRead"))) {
   methods::setAs(
      "character",
      "ShortRead",
      function(from) {
         if (length(from) == 1 && file.exists(from)) {
            from <- tryCatch(
               ShortRead::readFasta(from),
               error = function(e) {
                  ShortRead::readFastq(from)
               }
            )
         } else {
            ShortRead::ShortRead(
               sread = Biostrings::DNAStringSet(from, use.names = FALSE),
               id = Biostrings::BStringSet(names(from))
            )
         }
      }
   )
}

sreadq_to_qsDNAss <- function(from) {
   to = Biostrings::QualityScaledDNAStringSet(
      x = ShortRead::sread(from),
      quality = Biostrings::quality(from)
   )
   names(to) <- ShortRead::id(from)
   to
}


#' @importClassesFrom Biostrings QualityScaledXStringSet
#' @importClassesFrom ShortRead ShortReadQ

if (!methods::hasMethod("coerce", c("ShortReadQ", "QualityScaledXStringSet"))) {
   methods::setAs(
      "ShortReadQ",
      "QualityScaledXStringSet",
      sreadq_to_qsDNAss
   )
}

if (!methods::hasMethod("coerce", c("ShortReadQ", "QualityScaledDNAStringSet"))) {
   methods::setAs(
      "ShortReadQ",
      "QualityScaledDNAStringSet",
      sreadq_to_qsDNAss
   )
}

qsDNAss_to_sreadq <- function(seq) {
   ShortRead::ShortReadQ(
      sread = magrittr::set_names(methods::as(seq, "DNAStringSet"), NULL),
      quality = Biostrings::quality(seq),
      id = Biostrings::BStringSet(names(seq))
   )
}

#' @importClassesFrom Biostrings QualityScaledDNAStringSet

if (!methods::hasMethod("coerce", c("QualityScaledDNAStringSet", "ShortReadQ"))) {
   methods::setAs(
      "QualityScaledDNAStringSet",
      "ShortReadQ",
      qsDNAss_to_sreadq
   )
}


if (!methods::hasMethod("coerce", c("QualityScaledDNAStringSet", "ShortRead"))) {
   methods::setAs(
      "QualityScaledDNAStringSet",
      "ShortRead",
      qsDNAss_to_sreadq
   )
}