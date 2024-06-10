# TFM - Comparación de metodologías de agrupamiento en Transcriptómica Espacial

Datos, código e información del software empleado en el trabajo de fin de máster. Los documentos `.R` están pensados para ser ejecutados en segundo plano, pues su tiempo de ejecución puede superar las 24h en un hardware personal.

## Datos brutos

El conjunto de datos utilizado proviene del artículo [Maynard et al. 2021](https://www.nature.com/articles/s41593-020-00787-0).
Se puede descargar fácilmente en formato `SpatialExperiment` gracias a su paquete de `R`:

``` r
## Instala el paquete
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("spatialLIBD")

## Carga el paquete
library("spatialLIBD")

## Descarga los datos 
spe <- fetch_data(type = "spe")
```

## Preprocesamiento

El análisis previo a la ejecución de los métodos de agrupamiento sigue los pasos del libro online [Best Practices for Spatial Transcriptomics Analysis with Bioconductor](https://lmweber.org/BestPracticesST/). Está disponible en los documentos `qc.Rmd` y `Ploteros.Rmd`.

## Agrupamiento

La ejecución de los diferentes modelos de agrupamiento se llevó a cabo en el documento `clustering_cipf.R`. 

## Métricas de comparación

Las métricas de comparación se calcularon en el documento `tebanco.R`.

## Versión de R y paquetes empleados

``` r
sessionInfo()
```
R version 4.4.0 (2024-04-24 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 10 x64 (build 19045)

Matrix products: default


locale:
[1] LC_COLLATE=Spanish_Spain.utf8 
[2] LC_CTYPE=Spanish_Spain.utf8   
[3] LC_MONETARY=Spanish_Spain.utf8
[4] LC_NUMERIC=C                  
[5] LC_TIME=Spanish_Spain.utf8    

time zone: Europe/Madrid
tzcode source: internal

attached base packages:
[1] stats4    stats     graphics  grDevices utils    
[6] datasets  methods   base     

other attached packages:
 [1] ggrepel_0.9.5               dplyr_1.1.4                
 [3] ggspavis_1.10.0             reshape2_1.4.4             
 [5] patchwork_1.2.0             pheatmap_1.0.12            
 [7] scater_1.32.0               scran_1.32.0               
 [9] scuttle_1.14.0              cluster_2.1.6              
[11] ggplot2_3.5.1               Seurat_5.0.3               
[13] SeuratObject_5.0.1          sp_2.1-4                   
[15] spatialLIBD_1.16.0          SpatialExperiment_1.14.0   
[17] SingleCellExperiment_1.26.0 SummarizedExperiment_1.34.0
[19] Biobase_2.64.0              GenomicRanges_1.56.0       
[21] GenomeInfoDb_1.40.0         IRanges_2.38.0             
[23] S4Vectors_0.42.0            BiocGenerics_0.50.0        
[25] MatrixGenerics_1.16.0       matrixStats_1.3.0          

loaded via a namespace (and not attached):
  [1] spatstat.sparse_3.0-3     bitops_1.0-7             
  [3] httr_1.4.7                RColorBrewer_1.1-3       
  [5] doParallel_1.0.17         tools_4.4.0              
  [7] sctransform_0.4.1         utf8_1.2.4               
  [9] R6_2.5.1                  DT_0.33                  
 [11] uwot_0.2.2                lazyeval_0.2.2           
 [13] withr_3.0.0               gridExtra_2.3            
 [15] progressr_0.14.0          cli_3.6.2                
 [17] spatstat.explore_3.2-7    fastDummies_1.7.3        
 [19] sass_0.4.9                spatstat.data_3.0-4      
 [21] ggridges_0.5.6            pbapply_1.7-2            
 [23] Rsamtools_2.20.0          parallelly_1.37.1        
 [25] sessioninfo_1.2.2         attempt_0.3.1            
 [27] maps_3.4.2                limma_3.60.0             
 [29] rstudioapi_0.16.0         RSQLite_2.3.6            
 [31] generics_0.1.3            BiocIO_1.14.0            
 [33] spatstat.random_3.2-3     ica_1.0-3                
 [35] Matrix_1.7-0              ggbeeswarm_0.7.2         
 [37] fansi_1.0.6               abind_1.4-5              
 [39] lifecycle_1.0.4           yaml_2.3.8               
 [41] edgeR_4.2.0               SparseArray_1.4.0        
 [43] BiocFileCache_2.12.0      Rtsne_0.17               
 [45] paletteer_1.6.0           grid_4.4.0               
 [47] blob_1.2.4                dqrng_0.3.2              
 [49] promises_1.3.0            ExperimentHub_2.12.0     
 [51] crayon_1.5.2              miniUI_0.1.1.1           
 [53] lattice_0.22-6            beachmat_2.20.0          
 [55] cowplot_1.1.3             KEGGREST_1.44.0          
 [57] magick_2.8.3              metapod_1.12.0           
 [59] pillar_1.9.0              knitr_1.46               
 [61] rjson_0.2.21              future.apply_1.11.2      
 [63] codetools_0.2-20          leiden_0.4.3.1           
 [65] glue_1.7.0                data.table_1.15.4        
 [67] vctrs_0.6.5               png_0.1-8                
 [69] spam_2.10-0               gtable_0.3.5             
 [71] rematch2_2.1.2            cachem_1.0.8             
 [73] xfun_0.43                 S4Arrays_1.4.0           
 [75] mime_0.12                 ggside_0.3.1             
 [77] survival_3.5-8            iterators_1.0.14         
 [79] fields_15.2               bluster_1.14.0           
 [81] statmod_1.5.0             fitdistrplus_1.1-11      
 [83] ROCR_1.0-11               nlme_3.1-164             
 [85] bit64_4.0.5               filelock_1.0.3           
 [87] RcppAnnoy_0.0.22          bslib_0.7.0              
 [89] irlba_2.3.5.1             vipor_0.4.7              
 [91] KernSmooth_2.23-22        colorspace_2.1-0         
 [93] DBI_1.2.2                 tidyselect_1.2.1         
 [95] bit_4.0.5                 compiler_4.4.0           
 [97] curl_5.2.1                BiocNeighbors_1.22.0     
 [99] DelayedArray_0.30.0       plotly_4.10.4            
[101] rtracklayer_1.64.0        scales_1.3.0             
[103] lmtest_0.9-40             rappdirs_0.3.3           
[105] goftest_1.2-3             stringr_1.5.1            
[107] digest_0.6.35             spatstat.utils_3.0-4     
[109] rmarkdown_2.27            benchmarkmeData_1.0.4    
[111] XVector_0.44.0            htmltools_0.5.8.1        
[113] pkgconfig_2.0.3           sparseMatrixStats_1.16.0 
[115] dbplyr_2.5.0              fastmap_1.1.1            
[117] rlang_1.1.3               htmlwidgets_1.6.4        
[119] UCSC.utils_1.0.0          shiny_1.8.1.1            
[121] DelayedMatrixStats_1.26.0 jquerylib_0.1.4          
[123] zoo_1.8-12                jsonlite_1.8.8           
[125] BiocParallel_1.38.0       config_0.3.2             
[127] BiocSingular_1.20.0       RCurl_1.98-1.14          
[129] magrittr_2.0.3            GenomeInfoDbData_1.2.12  
[131] dotCall64_1.1-1           munsell_0.5.1            
[133] Rcpp_1.0.12               viridis_0.6.5            
[135] reticulate_1.36.1         stringi_1.8.3            
[137] zlibbioc_1.50.0           MASS_7.3-60.2            
[139] AnnotationHub_3.12.0      plyr_1.8.9               
[141] parallel_4.4.0            listenv_0.9.1            
[143] deldir_2.0-4              Biostrings_2.72.0        
[145] splines_4.4.0             tensor_1.5               
[147] locfit_1.5-9.9            igraph_2.0.3             
[149] spatstat.geom_3.2-9       RcppHNSW_0.6.0           
[151] ScaledMatrix_1.12.0       BiocVersion_3.19.1       
[153] XML_3.99-0.16.1           evaluate_0.23            
[155] golem_0.4.1               BiocManager_1.30.23      
[157] foreach_1.5.2             httpuv_1.6.15            
[159] polyclip_1.10-6           RANN_2.6.1               
[161] tidyr_1.3.1               purrr_1.0.2              
[163] future_1.33.2             benchmarkme_1.0.8        
[165] scattermore_1.2           rsvd_1.0.5               
[167] xtable_1.8-4              restfulr_0.0.15          
[169] RSpectra_0.16-1           later_1.3.2              
[171] viridisLite_0.4.2         tibble_3.2.1             
[173] memoise_2.0.1             beeswarm_0.4.0           
[175] AnnotationDbi_1.66.0      GenomicAlignments_1.40.0 
[177] shinyWidgets_0.8.6        globals_0.16.3  



