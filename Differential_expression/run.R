library(tidyverse)
library(ggrepel)

args=commandArgs(T)



d=read.table(as.character(args[1]),sep='\t',fill=T,quote="",header=T) %>% tibble
adjPvalCutoff=as.numeric(args[2])
plotTitle=as.character(args[3])
pal.sig=c("Significant"="#bababa","Not significant"="#404040", "Selected"="#ca0020")
geneList=c(	"TET1",
		"TET2",
		"TET3",
		"OGA",
		"OGT",
		"KMT2A",
		"KMT2B",
		"KMT2C",
		"KMT2D",
		"KMT2F",
		"KMT2G",
		"SETD1A",
		"SETD1B")

outprefix=as.character(args[4])

print("Checking gene symbols...")
for(gene in geneList){
	allGeneNames=d %>% pull(Gene.symbol)
	if(any(gene %in% allGeneNames)){
		print(paste0("Gene ", gene," found!"))
	}else{
		print(paste0("Gene ",gene," not found in the input table."))
	}
}

print("Analyzing for redundant proves...")
minPvalForGene=d %>% filter(Gene.symbol %in% geneList) %>% dplyr::select(ID,Gene.symbol,adj.P.Val) %>% group_by(Gene.symbol) %>% summarise(minPval=min(adj.P.Val)) %>% pull(minPval,Gene.symbol)

pdf(paste0(outprefix,".volcanoPlot.pdf"),height=5,width=4)
d4plot=d %>% mutate(label=ifelse(adj.P.Val<adjPvalCutoff,ifelse(Gene.symbol %in% geneList, ifelse(adj.P.Val==minPvalForGene[Gene.symbol],Gene.symbol,""),""),""), Significant=ifelse(adj.P.Val<adjPvalCutoff,"Significant","Not significant")) #%>% mutate(Significant=ifelse(Gene.symbol %in% geneList,"Selected",`Significant`)) 

d4plot.selected=d4plot %>% filter(label!="")
p=d4plot %>%
	ggplot(aes(x=logFC,y=-log(adj.P.Val),label=label,colour=Significant))+
	geom_point(size=0.5)+
	scale_color_manual(values=pal.sig)+
	geom_text(color="black",hjust = 0, nudge_x = 0.05)+
	#geom_text_repel(color="black",max.overlaps = Inf,force=1, force_pull=1,nudge_x=1,nudge_y=1,box.padding=0.5)+
	geom_point(data=d4plot.selected,colour="#ca0020",size=0.8)+
	theme_classic()+
	geom_hline(yintercept=-log(adjPvalCutoff),linetype="dashed",color="red")+
	geom_vline(xintercept=2,linetype="dashed",color="#67a9cf")+
	geom_vline(xintercept=-2,linetype="dashed",color="#67a9cf")+
	ggtitle(plotTitle)+
	theme(legend.title=element_blank(),legend.position="none")+
	xlab("Log2(Fold-change)")+
	ylab("-Log(Adj. P-Value)")
dev.off()
save.image(paste0(outprefix,".run.RData"))
