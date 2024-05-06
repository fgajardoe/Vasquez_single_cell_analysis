library(tidyverse)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)


# # Cargamos los datos
ogt.bw=import("../GSE39154/GSM957084_E14_Ogt.bigWig")
IgG.bw=import("../GSE39154/GSM957085_E14_IgG.bigWig")


# # Primer acercamiento al dato (Estamos asumiendo que el score refleja la magnitud/validez de los peaks)

# OGT
ogt.bw %>% as.data.frame %>% tibble %>% pull(score) %>% summary
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   6.00    6.00   11.00   12.58   17.00 1715.00


# IgG
IgG.bw %>% as.data.frame %>% tibble %>% pull(score) %>% summary
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   4.00    4.00    9.00   10.66   13.00 1355.00 


# # Creamos un ID para cada rango en el objeto GenomicRanges

IgG.bw=IgG.bw %>% as.data.frame %>% tibble %>% mutate(ID=paste0("IGG",1:NROW(.)))
ogt.bw=ogt.bw %>% as.data.frame %>% tibble %>% mutate(ID=paste0("OGT",1:NROW(.)))



# # Aplicar filtro por score
scoreCutoff=6
IgG.valid_ranges=IgG.bw %>% as.data.frame %>% tibble %>% filter(score>scoreCutoff) %>% pull(ID)
ogt.valid_ranges=ogt.bw %>% as.data.frame %>% tibble %>% filter(score>scoreCutoff) %>% pull(ID)




save.image("run.RData")
