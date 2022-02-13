
## inferorg: automatic inference and conversion of gene IDs 

`inferorg` is an `R` package to make it easy to find out the organism and gene ID format that a given set of IDs corresponds to. It is ideal for other packages that want to accept a set of gene IDs without requiring a particular format.

### Usage

To find out the organism and gene ID format, use the `inferorg` function:

```{r}
human_mhc_genes <- c("HLA-A", "HLA-B", "HLA-C")
inferorg(human_mhc_genes)
```
which outputs:
```{bash}
## $organism
## [1] "human"
## 
## $format
## [1] "symbol"
## 
## $confidence_organism
## [1] 1
## 
## $confidence_format
## [1] 1
```

The confidence scores provide an idea of certainty of the inferred organism and format: they are high only if the IDs map to that gene/format more selectively than any other gene/format, and if most of the IDs map to that gene/format. See the vignette for more details.

The `autoconvert` function will also convert between formats without having to specify the input format and organism:

```{r}
autoconvert(c(21354, 21355), to = "symbol")
```

outputs:

```{bash}
## [1] "Tap1" "Tap2"
```

### Installation

To install the development branch, run:

```{r}
devtools::install_github("camlab-bioml/inferorg")
```


