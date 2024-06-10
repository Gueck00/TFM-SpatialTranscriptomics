# TFM - Comparación de metodologías de agrupamiento en Transcriptómica Espacial

## Datos brutos

El conjunto de datos utilizado proviene del artículo [Maynard et al. 2021](https://www.nature.com/articles/s41593-020-00787-0).
Se puede descargar fácilmente gracias a su paquete de `R` con:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("spatialLIBD")

## Carga el paquete
library("spatialLIBD")

## Descarga los datos 
spe <- fetch_data(type = "spe")
```

