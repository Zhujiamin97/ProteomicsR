# ProteomicsR 

A tool for simplifying proteomics data analysis

## How to install?

``` r
if(!require(devtools)){
install.packages("devtools")
}

devtools::install_github("Zhujiamin97/ProteomicsR")
```

## 修饰项目鉴定深度

``` r
ProteomicsR::spectronaut_PTM_depth (rm_modify = c("C","M"),
                                    SiteProbability = 0.75)

```

## HELA质控

``` r
ProteomicsR::hela_QC_maxquant()
```
