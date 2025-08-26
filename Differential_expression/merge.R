load("Haferlach.run.RData")
p1=p
load("Homminga.run.RData")
p2=p
load("Cramer-morales.run.RData")
p3=p

library(grid)
library(gridExtra)

pdf("panel.pdf",height=4,width=12)
grid.arrange(p1,p2,p3,nrow=1)
dev.off()
