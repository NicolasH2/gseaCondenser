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

`condenseGsea` takes a data.frame where every row is a pathway and there must have a column with genes (each item in that row is a string where genes are separated by commas). Those could be all genes that were found or all genes that are in the pathway or just the leading edge. condenseGsea will return the same data.frame but with additional columns that give information about the parent-child relation between terms. Terms are related if they share enough genes with another pathway, namely if the percentage of shared genes crosses the `similarity` parameter. The smaller pathway (in terms of how many genes it has) will normally be the child of the bigger pathway.

``` r
library(gseaCondenser)

gsea <- gseaCondenser::gseaTest                                # data set built into this package for demonstration purposes
gsea0 <- condenseGsea(gsea, colname="genes", similarity=0.1, idcol="pathway")   # here we use a very low similarity threshold. 0.8-1 might be more appropriate in many cases
```

The data.frame now contains few additional columns:

- cID: arbitrary number to identify the pathway
- cKids: cID of the pathways that are considered children of this pathway
- cParents: cID of the pathways that are considered parents of this pathway
- cMotherID: cID of the pathway that is considered the best parent (usually highest overlap)
- cShared: Ratio of the number of overlapping genes to the number of genes in the pathway
- cEve: boolean value, stating if the term has no parents (except itself)
- cMother: name of the mother pathway. This is only given if the user defines idcol

Let's plot the data to get an impression of the grouping. In [ggplot2](https://ggplot2.tidyverse.org/) we use the cMother column for the coloring. We see that defense response, cellular oxidant detoxification and antimicrobial humoral response are grouped under defense response. Also, platelet degranulation became the mother term for platelet degranulation and innate immune response.

```r
library(ggplot2)
plotme <- function(x){
  ggplot(x, aes(NES, pathway, size=pval, color=cMother)) +
    geom_point() +
    theme_classic()
}

plotme(gsea0)
```
<img src="readme_files/gsea_standard.png"/>

We can define how many parents should exist (maximum).

```r
gsea1 <- condenseGsea(gsea, similarity=0.1, idcol="pathway", n_finalParents=6)
plotme(gsea1)
```

<img src="readme_files/gsea_nparents.png"/>

Alternatively, we can manually define the parent pathways. This means that we don't have to pathway a similarity threshold. In any case the pathways will get one of the defined parents, namely the one with the highest overlap.

```r
gsea2 <- condenseGsea(gsea, idcol="pathway", finalParents = c("neutrophil degranulation",
                                                              "proteolysis", 
                                                              "adenylate cyclase-activating G protein-coupled receptor signaling pathway"))
plotme(gsea2)
```

<img src="readme_files/gsea_specterms.png"/>

Be aware that defining `n_finalParents` or `finalParents` will force pathways to be children of pathways they might normally not have as parents. The similarity threshold will be disregarded (although `n_finalParents` uses it at first to determine suitable parents, so it is still important there).
In the next codechunk, we subset the data.frame, so that only terms with a similarity of 0.02 (which is already low) are in the plot. This throws out some terms that have almost 0 similarity with their parents.

```r
plotme(subset(gsea2, cShared>.02))
```

<img src="readme_files/gsea_specterms_filter.png"/>

# freqGsea

`freqGsea` takes a vector of pathway names, cuts them up into pieces and returns a frequency table of all words from the pathway names. This also includes word-combinations of words that directly follow one another in a single pathway name. By default some words will be removed ("of","via",etc.) but you can switch this off, or even add terms that you might want to omit. In this example the [wordcloud2](https://github.com/Lchiffon/wordcloud2) package is used to vizualize the result.

``` r
library(gseaCondenser)
library(wordcloud2)

myfreq <- freqGsea(names(gseaCondenser::mySetlist)) # data set built into this package for demonstration purposes
wordcloud2(myfreq)
```
<img src="readme_files/wordcloud.png"/>
