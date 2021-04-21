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
  genes <- strsplit(gsea[,colname], split=sep)

  #ratiomat will become matrix where each row is a set and each column is a set. The cells will contain the percentage of overlap between the sets.
  ratiomat <- sapply(genes, function(x) sapply(genes, function(y) length(intersect(x,y)) )) #first, the cells contain absolute numbers
  ratiomat <- as.data.frame(ratiomat)
  rownames(ratiomat) <- NULL
  colnames(ratiomat) <- seq(nrow(gsea))
  ratiomat <- apply(ratiomat, 2, function(x) x/max(x)) #ratio is intersect of the sets vs. set size (of the column set)
  ratiomat <- as.data.frame(ratiomat)

  #each time a similarity value reaches the threshold, the smaller set's ID is added to "eaten" and to the "condenseChildren" column of the bigger set
  eaten <- NA
  gsea$condenseChildren <- ""
  gsea$condenseParents <- ""
  for(i in seq(nrow(ratiomat))) { #the ratio is intersect/n (number of genes of the column set), So j
    for(j in seq(ncol(ratiomat))){
      #if the ratio of intersect and set size (of set j) is bigger than the ratio of intersect and set size (of set i), it means that set j is smaller. Thus, it will get eaten
      if(ratiomat[i,j]>similarity & ratiomat[i,j] > ratiomat[j,i]){ #by using >, sets do not eat themselves, and it is always the smaller set that gets eaten
        gsea$condenseChildren[i] <- paste0(gsea$condenseChildren[i],",",j) #the set that does not get eaten gets the ID of j for its Children column
        gsea$condenseParents[j] <- paste0(gsea$condenseParents[j],",",i) # the set that does get eaten gets the ID of i for its Parent column
        eaten <- c(eaten,j) # just a vector, gathering all IDs from sets that were eaten
      }
    }
  }
  gsea$condenseChildren <- gsub("^,","",gsea$condenseChildren)
  gsea$condenseParents <- gsub("^,","",gsea$condenseParents)

  parents <- strsplit(sapply(gsea$condenseParents, function(x) x), split=",")
  parents <- lapply(parents, as.numeric)

  #the best parent will be determined by which parent shares the most genes with the child
  parent <- lapply(seq_along(parents), function(j) {
    curparents <- parents[[j]]
    if(length(curparents)==0){NA}else{ #for each row in gsea, extract the parents (j is the child)
      overlaps <- sapply(curparents, function(i) ratiomat[j,i]) #extract overlap values for all parents (normalized by the set size of the child)
      parent <- curparents[which(overlaps==max(overlaps))] #choose the parent with the highest overlap value. This could lead to several parents

      overlaps2 <- ratiomat[ rep(j,length(parent)) , parent ] #extract overlap values for all parents (normalized by the set size of the parent)
      parent <- ifelse(length(parent)==1, parent, parent[which(overlaps2==max(overlaps2))][1] ) #choose the parent with the highest overlap value. Still, this could leave several parents, which is why we just pick the first one
      return(parent)
    }
  })
  parent <- ifelse(is.na(parent), gsea$condenseID, parent) #if there is no parent, the child becomes its own parent
  gsea$condenseParentID <- unlist(parent)
  gsea$condenseParentName <- gsea$pathway[gsea$condenseParentID]
  gsea$condenseSurvive <- ifelse(gsea$condenseID %in% eaten, FALSE, TRUE)

  return(gsea)
}
