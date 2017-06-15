#' Multivariate genome-wide association scan
#' 
#' The function imports GenABEL (gwaa.data class) or DatABEL (.fv*) data formats
#' and performs multivariate test for each genetic variant using multivariate
#' analysis of variance (MANOVA).
#' 
#' @param x An object created by \code{\link{MultiLoad}}.
#' @param trait.idx A vector giving the indices of traits to be analyzed.
#' @param ... not used.
#' 
#' @note Either \code{gwaa.data} (for GenABEL data format) or the combination of 
#' \code{phenofile} and \code{genofile} (for DatABEL data format) has to be provided.
#' If all are provided, only \code{phenofile} and \code{genofile} will be used. When using
#' DatABEL format input, individual IDs in \code{phenofile} and \code{genofile} have to match!
#' 
#' @return The function returns a data frame containing the multi-trait GWAS results, where the row names are
#' the variants names. The column names are: variant name (\code{Marker}), allele frequency (\code{Freq}),
#' the smallest sample size of the traits (\code{N}), effect on the phenotype score (\code{Beta.S}, see reference),
#' standard error (\code{SE}), p-value (\code{P}), and the rest the coefficients to construct the phenotype score
#' (see reference).
#' 
#' @author Xia Shen
#' 
#' @references 
#' Xia Shen, ..., Jim Wilson, Gordan Lauc, Yurii Aulchenko (2015).
#' Multi-omic-variate analysis identified novel loci associated with 
#' compound N-Glycosylation of human Immunoglobulin G. \emph{Submitted}.
#' 
#' @seealso 
#' \code{\link{MultiLoad}}
#' 
#' @examples 
#' 
#' ## loading example gwaa.data in GenABEL
#' require(GenABEL)
#' data(ge03d2ex.clean)
#' 
#' ## running multivariate GWAS for 3 traits: height, weight, bmi
#' loaded <- MultiLoad(gwaa.data = ge03d2ex.clean, trait.cols = c(5, 6, 8), 
#'                     covariate.cols = c(2, 3))
#' 
#' ## running the multivariate GWAS 
#' res <- Multivariate(loaded)
#' 
#' @aliases Multivariate, multivariate
#' @keywords multivariate
#' 
`Multivariate` <- function(x, trait.idx = NULL, ...) {
	if (class(x) != 'multi.loaded') {
		stop('wrong data type!')
	}
	if (!is.null(trait.idx)) {
		idx <- as.numeric(t(cbind(trait.idx*2 - 1, trait.idx*2)))
		x$gwa <- x$gwa[,idx]
		x$cor.pheno <- x$cor.pheno[trait.idx,trait.idx]
	}
	if (any(is.na(rowSums(x$gwa)))) {
		naidx <- which(is.na(rowSums(x$gwa)))
		x$gwa <- x$gwa[-naidx,]
		x$cvg <- x$cvg[-naidx]
		x$snp <- x$snp[-naidx]
		x$nsnp <- x$nsnp - length(naidx)
	}
	npheno <- ncol(x$cor.pheno)
    cat('initializing analysis ...')
	res <- matrix(NA, x$nsnp, 5)
	coef <- matrix(NA, x$nsnp, npheno)
	res <- data.frame(Marker = rownames(x$gwa), Beta.S = NA, SE = NA, P = NA)
    rownames(res) <- rownames(coef) <- x$snp
    cat(' OK\n')
	cat('multivariate GWA scan ...')
	k <- x$nsnp
	m <- nrow(x$cor.pheno)
	n <- rep(as.integer(x$nid), k)
	betamat <- x$gwa[,seq(1, 2*m, 2)]
	scan <- .Fortran('MultiSummaryLoopPrecise', k = as.integer(k), m = as.integer(m), nn = as.numeric(n), 
			varg = x$cvg, betamat = betamat, R = x$cor.pheno, invR = solve(x$cor.pheno), D = matrix(0, m, m),
			sdY = diag(sqrt(x$var.pheno)), invsdY = diag(1/sqrt(x$var.pheno)), sY = sqrt(x$var.pheno), 
			b = betamat, s = betamat, pil = numeric(k), coef = betamat, ss = betamat, PACKAGE = "MultiABEL")
	cat(' OK\n')
	res$Beta.S <- scan$pil
	scan$fstat <- res$Beta.S/m/(1 - res$Beta.S)*(n - m - 1)
	res$P <- pf(scan$fstat, m, n - m - 1, lower.tail = FALSE)
	res$SE <- res$Beta.S/sqrt(qchisq(res$P, 1, lower.tail = FALSE))
	coef <- scan$coef
	colnames(coef) <- paste('coef.', colnames(x$cor.pheno), sep  = '')
	return(data.frame(res, coef))
}
