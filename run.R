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




# # Aplicar filtro por cromosoma (xq el Mac no se la puede con mas)
whichChr="chr10"
IgG.bw=IgG.bw[seqnames(IgG.bw)==whichChr]
ogt.bw=ogt.bw[seqnames(ogt.bw)==whichChr]

# # Creamos un ID para cada rango en el objeto GenomicRanges -> candidato a ser borrado
IgG.ids=IgG.bw %>% as.data.frame %>% tibble %>% mutate(ID=paste0("IGG",1:NROW(.))) %>% pull(ID)
names(IgG.bw)=IgG.ids
IgG.bw$ID=IgG.ids

ogt.ids=ogt.bw %>% as.data.frame %>% tibble %>% mutate(ID=paste0("OGT",1:NROW(.))) %>% pull(ID)
names(ogt.bw)=ogt.ids
ogt.bw$ID=ogt.ids

# # Aplicar filtro por score
# IgG
scoreCutoff=10
IgG.tibble=IgG.bw %>% as.data.frame %>% tibble 
Ninicial=IgG.tibble %>% NROW
IgG.tibble=IgG.tibble %>% filter(score>scoreCutoff)
Nfinal=IgG.tibble %>% NROW
removedPerc=100-(Nfinal*100/Ninicial)
cat(paste0("Removed percentage (IgG): ",removedPerc," %"))
IgG.ids.overScoreCutoff=IgG.tibble %>% pull(ID)

gc()

# OGT
scoreCutoff=10
ogt.tibble=ogt.bw %>% as.data.frame %>% tibble 
Ninicial=ogt.tibble %>% NROW
ogt.tibble=ogt.tibble %>% filter(score>scoreCutoff)
Nfinal=ogt.tibble %>% NROW
removedPerc=100-(Nfinal*100/Ninicial)
cat(paste0("Removed percentage (OGT): ",removedPerc," %"))
ogt.ids.overScoreCutoff=ogt.tibble %>% pull(ID)
gc()

# Filter GenomicRanges by Score
IgG.bw=IgG.bw[IgG.ids.overScoreCutoff]
ogt.bw=ogt.bw[ogt.ids.overScoreCutoff]




save.image("run.RData")
