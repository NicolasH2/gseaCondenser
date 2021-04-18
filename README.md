# gseaCondenser

An [R](https://www.r-project.org) package that aims to give better overviews over GSEA.

Table of contents:

- [Installation](#Installation)
- [condenseGsea](#condenseGsea)
- [freqGsea](#freqGsea)

# Installation
Install the package from the git repository:
``` r
devtools::install_github("nicolash2/ggbrace")
```

# condenseGSEA

condenseGsea takes a data.frame that must have a column with genes. Those could be all genes that were found or all genes that are in the set or just the leading edge. condenseGsea will return the same data.frame but with 3 additional columns, most important of all the column condenseDropout, which signals if this set is redundant. Terms are redundant if they share enough genes with another set. The smaller set (in terms of how many genes it has) will always be eaten up by the bigger set and become redunant. Which sets have been consumed by a set is listed in the condenseChildren column. The numbers in it refer to the condenseID that the set was given (arbitrary numbering).

``` r
library(gseaCondenser)

gsea <- gseaCondenser::myGsea                 # data set built into this package for demonstration purposes
gsea <- condenseGsea(gsea, similarity=0.3)   # here we use a very low similarity threshold. 0.8-1 might be more appropriate in many cases
head(gsea)
```

The data.frame can now be filtered for those terms that are FALSE for the condenseDropout column. This will yield only non-redundant terms and thus you have cleaned up your GSEA. You can always go back and look into the condenseChildren column to see which terms are behind the non-redundant ones.

# freqGsea

freqGsea takes a vector of set names, cuts them up into pieces and returns a frequency table of all words from the set names. This also includes word-combinations of words that directly follow one another in a single set name. By default some words will be removed ("of","via",etc.) but you can switch this off, or even add terms that you might want not to be shown. Here we then use the [wordcloud2](https://github.com/Lchiffon/wordcloud2) package to vizualize the result.

``` r
library(gseaCondenser)
library(wordcloud2)

myfreq <- freqGsea(names(gseaCondenser::mySetlist))
wordcloud2::wordcloud2(myfreq)
```

