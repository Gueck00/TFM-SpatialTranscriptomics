# TFM - Comparación de metodologías de agrupamiento en Transcriptómica Espacial

Datos, código e información del software empleado en el trabajo de fin de máster. Los documentos `.R` están pensados para ser ejecutados en segundo plano, pues su tiempo de ejecución puede superar las 24h, en un hardware personal.

## Datos brutos

El conjunto de datos utilizado proviene del artículo [Maynard et al. 2021](https://www.nature.com/articles/s41593-020-00787-0).
Se puede descargar fácilmente gracias a su paquete de `R` con:

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

El análisis previo a la ejecución de los métodos de agrupamiento sigue los pasos del libro online [Best Practices for Spatial Transcriptomics Analysis with Bioconductor](https://lmweber.org/BestPracticesST/). 

## Agrupamiento



## Métricas de comparación




