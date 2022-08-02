
library(ggplot2)

# load original results
df <- read.csv(file = "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/fasthplus_resuts.csv", header = TRUE)

pdf(file = "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/plots/06_fasthplus/fasthplus_original.pdf")
ggplot(data = df, aes(x = k, y = X1.h, group = 1)) +
    geom_line() +
    geom_point()
dev.off()


# load data after removed bad spots
df <- read.csv(file = "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/06_fasthplus/fasthplus_results.csv", header = TRUE)

pdf(file = "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/plots/06_fasthplus/fasthplus.pdf")
ggplot(data = df, aes(x = k, y = X1.h, group = 1)) +
    geom_line() +
    geom_point()
dev.off()

# load data after removed whtie matter
df <- read.csv(file = "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/06_fasthplus/fasthplus_results_no_WM.csv", header = TRUE)

pdf(file = "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/plots/06_fasthplus/fasthplus_no_WM_080222.pdf")
ggplot(data = df, aes(x = k, y = X1.h, group = 1)) +
    geom_line() +
    geom_point() +
    ylab(expression(1 - H^{
        "+"
    })) +
    theme_bw(base_size = 20) +
    geom_vline(xintercept = 16, color = "red", size = 0.75) +
    annotate(geom = "text", label = "k = 16", x = 18.7, y = 0.624, color = "red", size = 5) +
    geom_vline(xintercept = 4, color = "gray", size = 0.75, linetype = "dotted") +
    annotate(geom = "text", label = "k = 4", x = 6.9, y = 0.624, color = "gray", size = 5)
dev.off()
