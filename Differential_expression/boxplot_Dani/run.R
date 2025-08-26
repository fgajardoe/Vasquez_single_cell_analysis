library(tidyverse)
library(limma)


# Frase tipo:
# Our results when comparing the distribution of intensities between T-ALL and healthy individuals suggests that, in general T-ALL individuals have higher expression of KMT2A and TET1. This is in coherence with previous reports in other types of cancer.


# functions	
generateBoxPlot_old=function(d,eq,gene,comparison){
	eq.samples=eq %>% filter(Leukemia_class %in% comparison)
	samples=eq.samples %>% pull(ID_REF)
	d.selectedCols=d %>% select(all_of(c("ID_REF",samples)))
	d.selectedGene=d.selectedCols %>% filter(ID_REF==gene) %>% t %>% as.data.frame %>% tibble
	d.selectedGene=d.selectedGene[-1,] %>% mutate(ID_REF=samples)
	#return(d.selectedGene)
	colnames(d.selectedGene)=c("Intensity","ID_REF")
	d.selectedGene$Intensity=as.numeric(d.selectedGene$Intensity)
	ret=eq.samples %>% left_join(d.selectedGene,by="ID_REF")
	return(ret)
}

generateBoxPlot=function(d,eq,gene){
	samples=colnames(d)[-1]
	eq.samples=eq %>% filter(ID_REF %in% samples)
	d.selectedGene=d %>% filter(ID_REF==gene) %>% t %>% as.data.frame %>% tibble
	d.selectedGene=d.selectedGene[-1,] %>% mutate(ID_REF=samples)
	colnames(d.selectedGene)=c("Intensity","ID_REF")
	d.selectedGene$Intensity=as.numeric(d.selectedGene$Intensity)
	ret=eq.samples %>% left_join(d.selectedGene,by="ID_REF")
	return(ret)
}

filterMatrix=function(d,eq,comparison){
	eq.samples=eq %>% filter(Leukemia_class %in% comparison)
	samples=eq.samples %>% pull(ID_REF)
	d.selectedCols=d %>% dplyr::select(all_of(c("ID_REF",samples)))
	return(d.selectedCols)
}

# Parameters
rawIntensities=T


# Main # # #
print("loading samples information...")
eq=read.table("GSE13204-GPL570_series_matrix.sampleMetadata.tsv",sep="\t",fill=T) %>% tibble
eq=eq %>% t %>% as.data.frame %>% tibble
colnames(eq)=c("ID_REF","Sample_type","Leukemia_class")
eq=eq[-1,]
print("done.")

if(rawIntensities==T){
	# loading
	print("loading matrix...")
	d=read.table("GSE13204-GPL570_series_matrix.noMetadata.tsv",sep="\t",fill=T,header=T) %>% tibble
	print("done.")

	# filtering columns
	print("filtering matrix...")
	comparison=c("leukemia class: Non-leukemia and healthy bone marrow","leukemia class: T-ALL")
	d=filterMatrix(d,eq,comparison)
	print("done.")

	# normalizing
	print("normalizing and correcting...")
	d.ids=d[,1]
	d=backgroundCorrect(d[,-1], method="normexp")
	d=normalizeBetweenArrays(as.data.frame(d))
	d.ready=cbind(d.ids,d) %>% tibble
	print("done.")

	# saving
	print("saving normalized intensities...")
	save(d.ready, file="normIntensities.RData")
	print("done.")

}else{
	print("Normalized intensities matrix found. Using this instead.")
	load("normIntensities.RData")
}


print("ready to plot boxplots...")


pdf("run.pdf")

gene="228906_at" #"TET1"
toPlot=generateBoxPlot(d.ready,eq,gene)
toPlot %>% ggplot(aes(x=Leukemia_class, y=Intensity, fill=Leukemia_class))+geom_boxplot()+ggtitle("TET1")

gene="1559856_s_at" #KMT2A
toPlot=generateBoxPlot(d.ready,eq,gene)
toPlot %>% ggplot(aes(x=Leukemia_class, y=Intensity, fill=Leukemia_class))+geom_boxplot()+ggtitle("KMT2A")

gene="207563_s_at" #OGT
toPlot=generateBoxPlot(d.ready,eq,gene)
toPlot %>% ggplot(aes(x=Leukemia_class, y=Intensity, fill=Leukemia_class))+geom_boxplot()+ggtitle("OGT")


dev.off()
