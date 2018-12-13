testthat::context("Testing read-fasta")
if (basename(getwd()) == 'testthat') setwd('../..')  # workspace is reset per file


# dir('data-raw/aliscore/', full.names = T) %>% file.remove()
# file.copy('data-raw/aligned-multiline.fasta', 'data-raw/aliscore/test.fasta');

# testthat::test_that("Testing aliscore", {
#     testthat::expect_identical(aliscore('-i data-raw/aliscore/test.fasta', F, F), 0L);
# });

# testthat::test_that("Testing alicut", {
#     testthat::expect_identical(alicut('data-raw/aliscore', '-s'), 0L);
# });

