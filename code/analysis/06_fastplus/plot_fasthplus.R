
library(ggplot2)

#load data
df <- read.csv(file = "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/fasthplus_resuts.csv", header = TRUE)

pdf(file = "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/plots/fasthplus.pdf")
ggplot(data=df, aes(x=k, y=X1.h, group=1)) +
  geom_line()+
  geom_point()
dev.off()