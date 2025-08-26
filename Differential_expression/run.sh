Rscript_cmd="docker run -it -v `pwd`:/home/user fgajardoe/r-custom:latest Rscript "
$Rscript_cmd run.R GSE13204.T-ALL_vs_Healthy.top.table.tsv 0.05 "Haferlach et al., 2010" Haferlach
$Rscript_cmd run.R Homminga_et_al_2011_GSE26713.top.table.tsv 0.05 "Homminga et al., 2011" Homminga
#$Rscript_cmd run.R Fernandes_et_al_2015_GSE66638.top.table.tsv 0.05 "Fernandes et al., 2015" Fernandes
$Rscript_cmd run.R Cramer-morales_et_al_2013_GSE48558.top.table.tsv 0.05 "Cramer-morales et al., 2013" Cramer-morales
#$Rscript_cmd run.R Borssén_et_al_2013_GSE41621.top.table.tsv 0.05 "Borssén et al., 2013" Borssén
#$Rscript_cmd run.R Anderson_et_al_2007_GSE7186.top.table.tsv 0.05 "Anderson et al., 2007" Anderson
$Rscript_cmd merge.R
