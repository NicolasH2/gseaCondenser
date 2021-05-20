# devtools::document()

#' add redundancy information to a GSEA table
#'
#' @param gsea data.frame, where each row is a gene set and there is at least one column with a string for each row. This string contains all genes for that set
#' @param colname string, naming the column in which the genes are listed
#' @param sep character. What genes in the column string? E.g. if a cell of the column would look like "HOXA9,HOXA3", you need to set sep=","
#' @param similarity number, what has to be the minimum gene overlap between two terms in order for one to be regarded as redundant?
#' @param n_finalParents integer, how many parents should there be among all pathways?
#' @param finalParents string vector, name the pathways that should be parents; works only if idcol is defined
#' @param idcol string, name of the column that contains the pathway names; necessary for finalParents; if provided adds an additional column for the parent name
#' @return data.frame similar to the input, but with 3 added columns: cID has a simple numeric ID for each row. cKids lists the IDs of all sets that are condensed under this set. cParents lists the possible parents for this set. cMotherID lists the sets that is the best possible parent among the parents. cMother is optionally added if idcol is provided. cEve states if the set has no parents. cShared gives the ratio of the overlap size (how many genes are shared between mother and child) versus the number of genes in the child term.
#' @export
#' @examples
#' library(gseaCondenser)
#'
#' gsea <- gseaCondenser::gseaTest
#' gsea <- condenseGsea(gsea, similarity=0.1, idcol="pathway")
#' head(gsea)
condenseGsea <- function(gsea, colname="genes", sep=",", similarity=0.9, n_finalParents=NULL, finalParents=NULL, idcol=NULL){
  if(!is.null(finalParents)){
    message("condenseGsea: Similarity cutoff will be ignored, because finalParents was defined.")
    similarity <- -1
  }
  gsea$cID <- seq(nrow(gsea))
  genes <- strsplit(gsea[,colname], split=sep)

  #ratiomat will become matrix where each row is a set and each column is a set. The cells will contain the percentage of overlap between the sets.
  ratiomat <- sapply(genes, function(x) sapply(genes, function(y) length(intersect(x,y)) )) #first, the cells contain absolute numbers
  ratiomat <- as.data.frame(ratiomat)
  rownames(ratiomat) <- NULL
  colnames(ratiomat) <- seq(nrow(gsea))
  ratiomat <- apply(ratiomat, 2, function(x) x/max(x)) #ratio is intersect of the sets vs. set size (of the column set)
  ratiomat <- as.data.frame(ratiomat)

  # define possible parents (only needed in case n_finalParents is defined)
  possibleParents <- seq(nrow(gsea))
  if(!is.null(finalParents)){
    if(is.null(idcol)) stop("please provide a column name for the column containing the pathway names")
    possibleParents <- which(gsea[,idcol] %in% finalParents)
  }

  #===================================================================================
  # define parents and children
  #===================================================================================
  # main function: decide parentage
  gsea <- parentage(gsea=gsea, ratiomat=ratiomat, similarity=similarity, possibleParents=possibleParents,
                    finalParents=finalParents, n_finalParents=n_finalParents)

  #============================================
  # decide parentage again (only with top parents and with no similarity cutoff) if n_finalParents was specified
  if(!is.null(n_finalParents)){
    possibleParents <- tail(sort(table(unlist(strsplit(gsea$cParents, split=",")))), n_finalParents) # possibleParents are the current top parents
    similarity <- -1 #similarity is lifted to allow the possibleParents to take in all children
    gsea <- parentage(gsea=gsea, ratiomat=ratiomat, similarity=similarity, possibleParents=possibleParents,
                      finalParents=finalParents, n_finalParents=n_finalParents) #re-do the parentage process, so only possible parents are chosen
  }

  #============================================
  parents <- strsplit(sapply(gsea$cParents, function(x) x), split=",")
  parents <- lapply(parents, as.numeric)

  #===================================================================================
  # define the mothers (best parent)
  #===================================================================================
  #the best parent will be determined by which parent shares the most genes with the child
  parentdata <- lapply(seq_along(parents), function(j) {
    curparents <- parents[[j]]
    if(length(curparents)==0){NA}else{ #for each row in gsea, extract the parents (j is the child)
      overlaps <- sapply(curparents, function(i) ratiomat[j,i]) #extract overlap values for all parents (normalized by the set size of the child)
      mother <- curparents[which(overlaps==max(overlaps))] #choose the parent with the highest overlap value. This could lead to several mothers

      overlaps2 <- ratiomat[ rep(j,length(mother)) , mother ] #extract overlap values for all parents (normalized by the set size of the parent)
      mother <- ifelse(length(mother)==1, mother, mother[which(overlaps2==max(overlaps2))][1] ) #choose the mother with the highest overlap value. Still, this could leave several mothers, so we just pick the first one
      return(c(mother,max(overlaps)))
    }
  })
  mother <- unlist(lapply(parentdata, function(x) x[1]))
  mother <- ifelse(is.na(mother), gsea$cID, mother) #if there is no parent, the child becomes its own parent

  #=======================
  #add columns
  gsea$cMotherID <- unlist(mother)
  gsea$cShared <- unlist(lapply(parentdata, function(x) x[2]))
  gsea$cShared <- ifelse(is.na(gsea$cShared), 1, gsea$cShared) # if there is no ratio, it is because the term is its own parent
  gsea$cEve <- ifelse(gsea$cMotherID == gsea$cID, TRUE, FALSE)

  if(!is.null(idcol)) gsea$cMother <- gsea[gsea$cMotherID,idcol]

  return(gsea)
}

# decides which terms are children and which are parents
parentage <- function(gsea, ratiomat, similarity, possibleParents, finalParents, n_finalParents){
  #each time a similarity value reaches the threshold, the smaller set's ID is added to "born" and to the "cKids" column of the bigger set
  born <- NA
  gsea$cKids <- ""
  gsea$cParents <- ""
  for(i in seq(nrow(ratiomat))) { #the ratio is intersect/n (number of genes of the column set), So j
    for(j in seq(ncol(ratiomat))){
      #if the ratio of intersect and set size (of set j) is bigger than the ratio of intersect and set size (of set i), it means that set j is smaller. Thus, it will get born
      if(ratiomat[i,j]>similarity & #the similarity threshold has to be passed
         i %in% possibleParents  &  # i has to be a listed
         (
           ratiomat[i,j] > ratiomat[j,i] | #the parent has to be bigger than the child (and not the two cannot be the same)
           !j %in% possibleParents |       #if j is not listed, than i can be a parent even if it is smaller than j
           (i==j & (!is.null(finalParents) | !is.null(n_finalParents))) #normally a term does not become its own parent, but if it is one of the listed finalParents, it will be
         )
      ){ #by using >, sets do not eat themselves, and it is always the smaller set that gets born
        gsea$cKids[i] <- paste0(gsea$cKids[i],",",j) #the set that does not get born gets the ID of j for its Children column
        gsea$cParents[j] <- paste0(gsea$cParents[j],",",i) # the set that does get born gets the ID of i for its Parent column
        born <- c(born,j) # just a vector, gathering all IDs from sets that were born
      }
    }
  }
  gsea$cKids <- gsub("^,","",gsea$cKids)
  gsea$cParents <- gsub("^,","",gsea$cParents)

  return(gsea)
}
