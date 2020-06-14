test_seq <- c(seq1 = "AAAAAAAAAAGGGGGGGGGGTTTTTTTTTTCCCCCCCCCC")
test_pos <- tibble::tibble(
   seq_id = rep("seq1", 4),
   region = LETTERS[1:4],
   start = c(1L, 11L, 21L, 31L),
   end = c(10L, 20L, 30L, 40L)
)
test_result <- c(seq1 = "TTTTTTTTTT")

test_ShortRead <- methods::as(test_seq, "ShortRead")
test_DNAStringSet <- Biostrings::DNAStringSet(test_seq)
test_RNAStringSet <- Biostrings::RNAStringSet(test_DNAStringSet)
test_fasta <- tempfile("test", fileext = ".fasta.gz")
Biostrings::writeXStringSet(test_DNAStringSet, test_fasta)
test_outfasta <- tempfile("test", fileext = ".fasta.gz")
extract_region(test_fasta, test_pos, region = "C", outfile = test_outfasta)

test_qual <- rep(10L, nchar(test_seq))
test_QualityScaledDNAStringSet <- Biostrings::QualityScaledDNAStringSet(
   test_DNAStringSet,
   Biostrings::PhredQuality(test_qual)
)
test_QualityScaledRNAStringSet <- Biostrings::QualityScaledRNAStringSet(
   test_RNAStringSet,
   Biostrings::PhredQuality(test_qual)
)
test_ShortReadQ <- methods::as(test_QualityScaledDNAStringSet, "ShortReadQ")
test_fastq <- tempfile("test", fileext = ".fastq.gz")
Biostrings::writeQualityScaledXStringSet(test_QualityScaledDNAStringSet, test_fastq)
test_outfastq <- tempfile("test", fileext = ".fastq.gz")
extract_region(test_fastq, test_pos, region = "C", outfile = test_outfastq)

test_that("extract_regions works for all formats", {
   expect_identical(
      extract_region(test_seq, test_pos, region = "C"),
      test_result
   )
   expect_equal(
      methods::as(
         extract_region(test_ShortRead, test_pos, region = "C"),
         "character"
      ),
      test_result
   )
   expect_equal(
      methods::as(
         extract_region(test_DNAStringSet, test_pos, region = "C"),
         "character"
      ),
      test_result
   )
   expect_equal(
      methods::as(
         extract_region(test_RNAStringSet, test_pos, region = "C"),
         "character"
      ),
      chartr("T", "U", test_result)
   )
   expect_equal(
      methods::as(
         extract_region(test_ShortReadQ, test_pos, region = "C"),
         "character"
      ),
      test_result
   )
   expect_equal(
      methods::as(
         extract_region(test_QualityScaledDNAStringSet, test_pos, region = "C"),
         "character"
      ),
      test_result
   )
   expect_equal(
      methods::as(
         extract_region(test_QualityScaledRNAStringSet, test_pos, region = "C"),
         "character"
      ),
      chartr("T", "U", test_result)
   )
   expect_equal(
      methods::as(
         Biostrings::readDNAStringSet(test_outfasta),
         "character"
      ),
      test_result
   )
})
