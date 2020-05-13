# copy of coerce code from LSUx, since it is not required for tzara

#' @importClassesFrom Biostrings DNAStringSet
#' @importClassesFrom ShortRead ShortRead

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



methods::setAs(
   "ShortRead",
   "DNAStringSet",
   function(from) {
      out <- ShortRead::sread(from)
      names(out) <- as.character(ShortRead::id(from))
      out
   }
)



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



methods::setAs(
   "ShortRead",
   "character",
   function(from) {
      out <- as.character(ShortRead::sread(from))
      names(out) <- as.character(ShortRead::id(from))
      out
   }
)


sreadq_to_qsDNAss <- function(from) {
   quality <- Biostrings::quality(from)
   if (methods::is(quality, "FastqQuality")) {
      quality <- methods::as(quality, "PhredQuality")
   } else if (methods::is(quality, "SFastqQuality")) {
      quality <- methods::as(quality, "SolexaQuality")
   }
   to = Biostrings::QualityScaledDNAStringSet(
      x = ShortRead::sread(from),
      quality = quality
   )
   names(to) <- ShortRead::id(from)
   to
}


#' @importClassesFrom Biostrings QualityScaledXStringSet
#' @importClassesFrom ShortRead ShortReadQ

methods::setAs(
   "ShortReadQ",
   "QualityScaledXStringSet",
   sreadq_to_qsDNAss
)


methods::setAs(
   "ShortReadQ",
   "QualityScaledDNAStringSet",
   sreadq_to_qsDNAss
)


qsDNAss_to_sreadq <- function(from) {
   ShortRead::ShortReadQ(
      sread = magrittr::set_names(methods::as(from, "DNAStringSet"), NULL),
      quality = Biostrings::quality(from),
      id = Biostrings::BStringSet(names(from))
   )
}

#' @importClassesFrom Biostrings QualityScaledDNAStringSet


methods::setAs(
   "QualityScaledDNAStringSet",
   "ShortReadQ",
   qsDNAss_to_sreadq
)




methods::setAs(
   "QualityScaledDNAStringSet",
   "ShortRead",
   qsDNAss_to_sreadq
)