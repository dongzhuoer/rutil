

#' @title Read sequences from a fasta file.
#'
#' @description `read_fasta()` reads all sequences contained in a fasta file
#'   into a tibble (see **Value**). See _FASTA format_ for the requirements of
#'   input file.
#'
#' @details Previously, I implemented fasta as a named character (in
#'   [bioinfor::read_fasta()]). But I found it more and more inconvenient as I
#'   used it. Many times when you manipulate the character, the name got lost.
#'   Finally, I decided to reimplement it as a tibble. That why this package
#'   come into being.
#'
#'   Actually, `file` can be anything as long as it's recognized by
#'   [readr::read_lines()].
#'
#'   If `per_line` isn't `TRUE`, `read_fasta()` will check it by detecting '>'
#'   and count line number. This (implemented by [stringr::str_detect()]) takes
#'   almost the same time as reading file content (implemented by
#'   [readr::read_lines()]). So specify `TRUE` if you are sure to save time. But
#'   this checking doesn't waste any time when sequences indeed span multiple
#'   lines. Additionly, a file contains odd number of lines will cause
#'   `pre_line=TURE` to be ignored
#'
#'   For `unalign`, when deal with large aligned file with linebreak containing
#'   many '-'s, implement _unalign_ ([stringr::str_replace_all()] in
#'   [read_fasta()]) is about 15% faster than do that externally, i.e. after
#'   [read_fasta()] returns.
#'
#' @param file string. Path to the fasta file.
#' @param per_line logical scalar. Whether sequences keep in one line or might
#'   span multiple lines. Specify `FALSE` if not sure, refer to **Details**.
#' @param unalign logical scalar. Whether unalign aligned sequences, i.e. remove
#'   '-', see **Details**.
#'
#' @return tibble
#'
#' - name: character. sequence header (without '>')
#'
#' - seq:  character. sequence itself, I named it `seq` instead of `sequence` to
#' save 62.5% of time (and everyone using R should take it for granted that seq
#' means sequence)
#'
#'
#' @examples
#' system.file('extdata', 'example.fasta', package = 'biozhuoer') %>% read_fasta();
#'
#' system.file('extdata', 'aligned-multiline.fasta', package = 'biozhuoer') %>% read_fasta();
#'
#' system.file('extdata', 'aligned-multiline.fasta', package = 'biozhuoer') %>% read_fasta(unalign = TRUE);
#'
#'
#' # crazy examples
#' read_fasta('>na>me\nATCG');
#'
#' read_fasta('>name\nAT>CG');
#'
#' read_fasta('>\nATCG');
#'
#' read_fasta('>name\n')
#'
#' read_fasta(paste0(c('>', rep('x', 10000), '\nATCG'), collapse = ''));
#'
#' tempfile() %T>% readr::write_file('>name', .)  %T>% readr::read_file() %>% read_fasta;
#'
#'
#' @section FASTA format: The minimum requirement is that it can only contain
#'   sequences (can use any character except for '>', at least at the begining)
#'   and their headers (begin with '>', one line). Comments are not supported.
#'
#'   Basically, `c('>na', 'me', 'ATCG')` (name span two lines) and `c('>name',
#'   '>ATCG')` (sequence begins with '>') are most common (and maybe the only)
#'   errors. Expect that, almost anything is acceptable (though your file may
#'   not be recongnized by other program):
#'
#'   1. name contains '>', see crazy example 1.
#'
#'   2. sequence contains '>' (can't be at the begining), see crazy example 2.
#'
#'   3. empty name, see crazy example 3.
#'
#'   4. empty sequence, see crazy example 4. But it must occupy an empty line if
#'   you want to scape check for `pre_line` (refer to details).
#'
#'   5. very long name, see crazy example 5. It can contain thousands of
#'   thousands characters, so long as you have enough memory (and disk if you
#'   are really mad).
#'
#'   6. don't you think the above is enough? We can't make friends anymore.
#'   (Never tell me you have tried crazy example 6)
#'
#'   In most situations, you won't have a file like that, unless you
#'   deliberately create one.
#'
#' @section Test files:
#'
#'   1. test1.fasta big file with linebreak.
#'
#'   1. test12.fasta same file without.
#'
#'   1. test2.fasta big aligned file with linebreak.
#'
#' @export
#'
read_fasta <- function(file, per_line = FALSE, unalign = FALSE) {
	#file  = 'data-raw/test.fasta';
	content <- readr::read_lines(file);
	N <- length(content);

	if (N %% 2 == 1) per_line = FALSE;   #per_line must be FALSE in that case even user specify TRUE

	if (!per_line) {
		header_line <- stringr::str_detect(content, '^>');
		n <- sum(header_line);
		if (n*2 == N) per_line = TRUE;
	}

	if (per_line) {
		seq  <- content[seq(2, N, 2)];
		name <- content[seq(1, N, 2)];
	} else {
		seq  <- character(n);
		name <- content[header_line];

		content[header_line] = '>';
		if (unalign) content %<>% stringr::str_replace_all('-', '');
		content %<>% paste0(collapse = '');

		seq = stringr::str_split(content, '>')[[1]][-1]; #the first line is empty
	}

	name %<>% stringr::str_replace('^>', '');
	tibble::data_frame(name, seq);
}

#Rprof(interval = 0.0002);read_fasta('data-raw/test12.fasta', F);Rprof(NULL);summaryRprof()
#Rprof(interval = 0.0002);read_fasta('data-raw/test12.fasta', T);Rprof(NULL);summaryRprof()
#Rprof(interval = 0.0002);read_fasta('data-raw/test2.fasta', F, T);Rprof(NULL);summaryRprof()
#Rprof(interval = 0.0002);invisible(stringr::str_replace_all(read_fasta('data-raw/test2.fasta', F, F), '-', ''));Rprof(NULL);summaryRprof()



#' @title Write sequences to a fasta file.
#'
#' @description overwrite
#'
#' @param fasta tibble. see **Value** of [read_fasta()].
#' @param path string. Path to the fasta file which you want write to.
#' @param width logical scalar. NotYetUsed.
#'
#' @return the input invisibly
#'
#' @examples
#' {
#'     input_file  <- system.file('extdata', 'example.fasta', package = 'biozhuoer')
#'     output_file <- tempfile();
#'     write_fasta(read_fasta(input_file), output_file)
#' }
#'
#' @export

write_fasta <- function(fasta, path, width) {
	if (!missing(width)) .NotYetUsed(width);

	N <- nrow(fasta) * 2;
	content <- character(N);

	content[seq(1, N, 2)] = paste0('>', fasta$name);
	content[seq(2, N, 2)] = fasta$seq;

	readr::write_lines(content, path);

	invisible(fasta)
}
