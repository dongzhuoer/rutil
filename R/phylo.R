#' @title wrapper to run Aliscore
#'
#' @details aliscore will produce several files under the directory containg input files
#'
#' @param args character. options for `Aliscore.0x.y.pl`
#' @param ... other parameters passed on to [system2()]
#'
#' @return result of [system2()]
#' @export
#'
#' @examples
#' \donotrun {
#'     aliscore('-i data-raw/aliscore/test.fasta')
#' }
aliscore <- function(args, ...) {
	args %<>% c('-I', pkg_file('lib'), pkg_file('exec/Aliscore.02.2.pl'), .)
	system2('perl', args, ...);
}


#' @title wrapper to run Alicut
#'
#' @param dir string. the directory your input files lie. working directory will
#'   be temporarily set to it
#' @param ... other parameters passed on to [system2()], including `args` and so
#'   on.
#'
#' @return result of [system2()]
#' @export
#'
#' @examples
#' \donotrun {
#'     aliscore('-i data-raw/aliscore/test.fasta');
#'     alicut('data-raw/aliscore', args = '-s')
#' }
alicut <- function(dir, ...) {
	old_wd <- getwd();
	on.exit(setwd(old_wd));

	setwd(dir);
	system2(pkg_file('exec/ALICUT_V2.31.pl'), ...);
}




