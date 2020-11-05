test_st <-
    list(
        region1 = tibble::tribble(
            ~seq_id, ~dada_seq,
            "seq1",  "seq1region1;",
            "seq2",  "seq2region1;",
            "seq4",  "seq4region1;"
        ),
        region2 = tibble::tribble(
            ~seq_id, ~dada_seq,
            "seq1",  "seq1region2;",
            "seq2",  "seq2region2;",
            "seq3",  "seq3region2;"
        )
    )

outputlist1 <-
    list(
        full = c("region1", "region2")
    )

outputlist2 <-
    list(
        region1 = "region1",
        region2 = "region2"
    )

expected_val <- tibble::tribble(
    ~seq_id, ~region1,       ~region2,
    "seq1",  "seq1region1;", "seq1region2;",
    "seq2",  "seq2region1;", "seq2region2;",
    "seq4",  "seq4region1;", NA_character_,
    "seq3",  NA_character_,  "seq3region2;"
)

test_that("assemble_region_table works", {
    expect_equal(
        tzara:::assemble_region_table(test_st, output = outputlist1),
        expected_val
    )
    expect_equal(
        tzara:::assemble_region_table(test_st, output = outputlist2),
        expected_val
    )
})
