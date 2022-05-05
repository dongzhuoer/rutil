testthat::context("Testing read-fasta")
setwd(here::here(''))  # workspace is reset per file



testthat::test_that("Testing read_fasta", {
    testthat::expect_true(identical(
        read_fasta('>name\nATCG'),
        tibble::tibble(name = 'name', seq = 'ATCG')
    ));
    
	testthat::expect_true(identical(
        read_fasta('>name1\nATCG\n>name2\nAGTC'),
        tibble::tibble(name = c('name1', 'name2'), seq = c('ATCG', 'AGTC'))
    ))
});
