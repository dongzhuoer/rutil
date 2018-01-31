# context("Testing phylo")

# setwd('../..');

# dir('data-raw/aliscore/', full.names = T) %>% file.remove()
# file.copy('data-raw/aligned-multiline.fasta', 'data-raw/aliscore/test.fasta');

# test_that("Testing aliscore", {
#     expect_identical(aliscore('-i data-raw/aliscore/test.fasta', F, F), 0L);
# });

# test_that("Testing alicut", {
#     expect_identical(alicut('data-raw/aliscore', '-s'), 0L);
# });

