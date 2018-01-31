context("Testing read-fasta")

#setwd('../..');

test_that("Testing read_fasta", {
    expect_true(identical(
        read_fasta('>name\nATCG'),
        tibble::data_frame(name = 'name', seq = 'ATCG')
    ));
	expect_true(identical(
        read_fasta('>name1\nATCG\n>name2\nAGTC'),
        tibble::data_frame(name = c('name1', 'name2'), seq = c('ATCG', 'AGTC'))
    ))
});


