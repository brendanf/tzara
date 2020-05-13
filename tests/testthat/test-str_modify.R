test_that("str_modify works when no modification is given", {
   expect_identical(
      tzara:::str_modify("test string"),
      "test string"
   )
   expect_identical(
      tzara:::str_modify("test string", regex = NULL),

      "test string"
   )
   expect_identical(
      tzara:::str_modify("test string", replace = NULL),
      "test string"
   )
   expect_identical(
      tzara:::str_modify("test string", regex = NULL, replace = NULL),
      "test string"
   )
})

test_that("str_modify works in extract mode", {
   expect_identical(
      tzara:::str_modify("test string", regex = "test"),
      "test"
   )
   expect_identical(
      tzara:::str_modify("test string", regex = "test", replace = NULL),
      "test"
   )
   expect_identical(
      tzara:::str_modify("test string", regex = "[^ ]+", replace = NULL),
      "test"
   )
   expect_identical(
      tzara:::str_modify("test string", regex = "missing", replace = NULL),
      NA_character_
   )
})

test_that("str_modify works in replace mode", {
   expect_identical(
      tzara:::str_modify("test string", regex = "test", replace = "result"),
      "result string"
   )
   expect_identical(
      tzara:::str_modify("test string", regex = "([^ ]+)", replace = "\\1"),
      "test string"
   )
   expect_identical(
      tzara:::str_modify("test string", regex = "missing", replace = "result"),
      "test string"
   )
})