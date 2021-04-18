# devtools::document()

#' frequence table for GSEA terms
#'
#' @param pathways character vector. Upper or lower case makes no matter, since everything will be converted to uppercase anyway.
#' @param sep charachter. What separates words in a set name? E.g. if terms look like "KEGG_CELL_CYCLE", you need to set sep="_"
#' @param bias numeric vector of the same length as the pathways vector The pathways will be counted x times (x = bias). The numbers will be rounded to integers.
#' @param clean boolean, should common words and wordcombinations including these common words be removed? These common words include e.g. "the","of", etc. but also "GO","KEGG","REACTOME" and "PID"
#' @param cleanII boolean, should unspecific words (but not their combination with others) be removed? E.g. "cellular","activation","alpha" etc. Combinations like il_alpha will not be removed.
#' @param removeIncluding character vector. Listed words and combinations that include them will be removed.
#' @param removeThis character vector. Listed words will be removed, but not combinations that include them.
#' @return Data.frame containing one column "word" and one column "freq", showing how many times one word (or word-combination) came up.
#' @export
#' @examples
#' library(gseaCondenser)
#' library(wordcloud2)
#'
#' myfreq <- freqGsea(names(gseaCondenser::mySetlist))
#' wordcloud2::wordcloud2(myfreq)
freqGsea <- function(pathways, sep="_", bias=NULL, clean=T, cleanII=T, removeIncluding="", removeThis=""){

  if(!is.null(bias)) pathways <- rep(pathways, round(bias))

  allcombis <- function(pathway){
    words <- unlist(strsplit(pathway, sep))

    wordrange <- seq_along(words)
    combis <- lapply(wordrange, function(y){ #for each word y ...
      sapply(wordrange, function(x) if(y>=x) paste(words[x:y], collapse="-")) # ...get a sequence of words from a word (x) in the term to the word (y)
    })
    combis <- unlist(combis)
    return(combis)
  }

  allterms <- unlist(lapply(pathways, allcombis)) #takes longest
  allterms <- as.data.frame(table(allterms))
  colnames(allterms) <- c("word","freq")

  # remove certain terms
  if(clean){
    discard <- c("the","to","from","with","an","and","or","of","off","on","by","via","in","kegg","go","reactome","pid")
    discard <- c(discard,removeIncluding)
    discard <- c(paste0("^",discard,"-"),paste0("-",discard,"-"),paste0("-",discard,"$"),paste0("^",discard,"$"))
  }
  discard2 <- ""
  if(cleanII) discard2 <- c("a","activation","mediated","events","event","cell","cellular","activity","regulation","positive","negative","process","alpha","beta","gamma","delta",".","ii","iii","iv","v","vi","vii","viii","ix")
  discard2 <- c(discard2, removeThis)
  discard2 <- c(paste0("^",discard2,"$"))

  discard <- c(discard, toupper(discard), discard2, toupper(discard2))
  discard <- paste(discard, collapse="|")

  allterms <- allterms[!grepl(discard, allterms$word),]
  allterms <- allterms[order(allterms$freq, decreasing=T),]

  #output
  return(allterms)

}

