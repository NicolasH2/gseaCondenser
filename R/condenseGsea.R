# devtools::document()

#' add redundancy information to a GSEA table
#'
#' @param gsea data.frame, where each row is a gene set and there is at least one column with a string for each row. This string contains all genes for that set
#' @param colname string, naming the column in which the genes are listed
#' @param sep charachter. What genes in the column string? E.g. if a cell of the column would look like "HOXA9,HOXA3", you need to set sep=","
#' @param similarity number, what has to be the minimum gene overlap between two terms in order for one to be regarded as redundant?
#' @return data.frame similar to the input, but with 3 added columns: condenseID has a simple numeric ID for each row. condenseChildren lists the IDs of all sets that were eaten by this set. condenseDropout states whether or not this set was eaten itself.
#' @export
#' @examples
#' library(gseaCondenser)
#'
#' gsea <- gseaCondenser::myGsea
#' gsea <- condenseGsea(gsea, similarity=0.3)
#' head(gsea)
condenseGsea <- function(gsea, colname="genes", sep=",", similarity=0.9){
  gsea$condenseID <- seq(nrow(gsea))
  gsea$condenseChildren <- ""
  genes <- strsplit(gsea[,colname], split=sep)

  #ratiomat will become matrix where each row is a set and each column is a set. The cells will contain the percentage of overlap between the sets.
  ratiomat <- sapply(genes, function(x) sapply(genes, function(y) length(intersect(x,y)) ))
  ratiomat <- as.data.frame(ratiomat)
  rownames(ratiomat) <- NULL
  colnames(ratiomat) <- seq(nrow(gsea))
  ratiomat <- apply(ratiomat, 2, function(x) x/max(x))

  #each time a similarity value reaches the threshold, the smaller set's ID is added to "eaten" and to the "condenseChildren" column of the bigger set
  eaten <- NA
  for(i in seq(nrow(ratiomat))) {
    for(j in seq(ncol(ratiomat))){
      if(ratiomat[i,j]>similarity & ratiomat[i,j] > ratiomat[j,i]){ #by using >, sets do not eat themselves, and it is always the smaller set that gets eaten
        gsea$condenseChildren[i] <- paste0(gsea$condenseChildren[i],",",j)
        eaten <- c(eaten,j) # just a vector, gathering all IDs from sets that were eaten
      }
    }
  }

  gsea$condenseChildren <- gsub("^,","",gsea$condenseChildren)
  gsea$condenseDropout <- ifelse(gsea$condenseID %in% eaten, TRUE, FALSE)
  return(gsea)
}
