#' Colon Cancer Dataset
#' 
#' The dataset describes a colon cancer study (\cite{Alon et al., 1999}) in
#' which gene expression levels were measured on 40 normal tissues and on 22
#' tumor colon tissues. Note that a logarithmic (base 10) transformation has
#' been applied to the gene expression levels.
#'
#' @name colon
#' @docType data
#' @format A data frame in which the rows correspond to 2000 genes and the
#' columns to 62 tissues. The first 40 columns belong to the normal tissue
#' group while the last 22 columns to the tumor colon tissue group.
#' @references Alon, U., Barkai, N., Notterman, D.A., Gish, K., Ybarra, S.,
#' Mack, D. and Levine, A.J. (1999) Broad patterns of gene expression revealed
#' by clustering analysis of tumor and normal colon tissues probed by
#' oligonucleotide arrays. \emph{Proceedings of the National Academy of
#' Sciences of the United States of America} \bold{96}, 6745--6750.
#' @source
#' \href{http://genomics-pubs.princeton.edu/oncology/affydata}{http://genomics-pubs.princeton.edu/oncology/affydata}
#' [Last Assessed: 2016-05-21]
#' @keywords datasets
#' @examples
#' data(colon)
#' summary(colon)
NULL
