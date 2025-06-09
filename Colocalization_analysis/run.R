library(tidyverse)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)

# Functions # # #

`%notin%` <<- Negate(`%in%`)

preprocessing=function(GRanges,
		       chr="chr10",
		       scoreCutoff=10,
		       prefix="ID"){



	# # Aplicar filtro por cromosoma (xq el Mac no se la puede con mas)
	whichChr=chr
	GRanges=GRanges[seqnames(GRanges)==whichChr]

	# # Creamos un ID para cada rango en el objeto GenomicRanges -> candidato a ser borrado
	GRanges.ids=GRanges %>% as.data.frame %>% tibble %>% mutate(ID=paste0(prefix,1:NROW(.))) %>% pull(ID)
	names(GRanges)=GRanges.ids
	GRanges$ID=GRanges.ids

	# # Aplicar filtro por score
	GRanges.tibble=GRanges %>% as.data.frame %>% tibble 
	Ninicial=GRanges.tibble %>% NROW
	GRanges.tibble=GRanges.tibble %>% filter(score>scoreCutoff)
	Nfinal=GRanges.tibble %>% NROW
	removedPerc=100-(Nfinal*100/Ninicial)
	cat(paste0("Removed percentage (",prefix,"): ",removedPerc," %"))
	GRanges.ids.overScoreCutoff=GRanges.tibble %>% pull(ID)
	gc()

	# Filter GenomicRanges by Score
	GRanges.len.i=GRanges %>% length
	GRanges=GRanges[GRanges.ids.overScoreCutoff]
	GRanges.len.f=GRanges %>% length
	GRanges.percentFiltered=(GRanges.len.i-GRanges.len.f)*100/GRanges.len.i

	# merge redundant ranges
	GRanges=GenomicRanges::reduce(GRanges)

	# nuevos IDs para los rangos en el objeto reducido
	GRanges.ids=GRanges %>% as.data.frame %>% tibble %>% mutate(ID=paste0(prefix,"R",1:NROW(.))) %>% pull(ID)
	names(GRanges)=GRanges.ids
	GRanges$ID=GRanges.ids


	return(GRanges)
}

removeNegativeControlOverlaps=function(GRanges, NegativeControlGRanges){
	# calculate overlapping
	o=findOverlapPairs(GRanges, NegativeControlGRanges)
	GRanges.unwantedIDs=o@first %>% data.frame %>% tibble %>% pull(ID)

	# exclude overlapping ranges. Keep only (OGT) ranges that do not overlap IgG.
	GRanges=GRanges[GRanges$ID %notin% GRanges.unwantedIDs]
	return(GRanges)
}



# # MAIN # #


# # Cargamos los datos
ogt.bw=import("../GSE39154/GSM957084_E14_Ogt.bigWig")
IgG.bw=import("../GSE39154/GSM957085_E14_IgG.bigWig")
ogt.bw=preprocessing(ogt.bw,prefix="OGT",scoreCutoff=10)
IgG.bw=preprocessing(IgG.bw,prefix="IGG",scoreCutoff=10)
ogt.bw=removeNegativeControlOverlaps(ogt.bw, IgG.bw)

rm(IgG.bw)
gc()

shSCR_TET1_C=rtracklayer::import("../GSM611192_Tet1_C_shSCR.bowtie.bed.gz")
IgG_TET1=rtracklayer::import("../GSM611198_IgG_ChIP.bowtie.bed.gz")
shSCR_TET1_C=preprocessing(shSCR_TET1_C,prefix="TET1_C",scoreCutoff=10)
IgG_TET1=preprocessing(IgG_TET1,prefix="IGG_TET1",scoreCutoff=10)
shSCR_TET1_C=removeNegativeControlOverlaps(shSCR_TET1_C, IgG_TET1)

rm(IgG_TET1)
gc()

library(TxDb.Mmusculus.UCSC.mm9.knownGene)
txdb=TxDb.Mmusculus.UCSC.mm9.knownGene

promoters=txdb %>% promoters

#o.counts=countOverlaps(ogt.bw,promoters)
o.counts=countOverlaps(promoters, ogt.bw)

save.image("run.RData")
