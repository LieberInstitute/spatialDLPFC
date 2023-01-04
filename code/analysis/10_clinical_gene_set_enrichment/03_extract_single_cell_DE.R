
library("xlsx")
library("here")
library("tidyverse")

dir_rdata <- here::here(
  "processed-data",
  "rdata",
  "spe",
  "10_clinical_gene_set_enrichment"
)

list.files(dir_rdata)

#### Nagy ####
# -Nagy et al.: Table 6 is DE by clusters (bimodal); Tables 8-31 are DEGs in each cluster; Table 32 is DE results with FDR correction across all clusters.  I think here we want Tables 8-31?
fn_nagy <- here(dir_rdata, "Nagy_MDD", "Nagy_41593_2020_621_MOESM3_ESM.xlsx")

# n <- 15
# # map(5:28, function(n){
# map(8:28, function(n){
#   data <- readxl::read_excel(fn_nagy, sheet = n)
#   header <- colnames(data_nagy)
#   ct <- gsub(".*results for (\\w.*) \\(.*","\\1",header[[1]])
#   
#   colnames(data) <- data_batiuk[2,]
#   data <- tail(data, -2) |> mutate(cellType = ct)
#   
#   return(data)
# })

## I think this has all the info we need
data_nagy <- readxl::read_excel(fn_nagy, sheet = "Supplementary Table 32", skip = 2)

# Gene     Estimate `Std. Error`    df `t value`    `Pr(>|t|)` `Fold Change`      FDR Cluster
# <chr>       <dbl>        <dbl> <dbl>     <dbl>         <dbl>         <dbl>    <dbl> <chr>  
# 1 TMSB4X    -0.402        0.0658  646.     -6.10 0.00000000178         0.669 0.000154 Astro3 
# 2 ACAA1     -0.113        0.0188  963.     -6.01 0.00000000266         0.893 0.000154 Ex7    
# 3 RPS11     -0.0800       0.0143 3368.     -5.58 0.0000000252          0.923 0.000977 Ex7    
# 4 KAZN       0.485        0.0958  263.      5.06 0.000000780           1.62  0.0211   OPC2   
# 5 PRAF2     -0.0616       0.0123  383.     -4.99 0.000000907           0.940 0.0211   Ex3    

data_nagy |> 
  filter(FDR < 0.05)  |>
  count(Cluster)

# Cluster     n
# <chr>   <int>
# 1 Astro3      2
# 2 Endo        2
# 3 Ex2         1
# 4 Ex3         2
# 5 Ex6         1
# 6 Ex7         3
# 7 Ex8         1
# 8 In1         1
# 9 In2         4
# 10 In8         1
# 11 OPC2        3
 
#### Velmeshev ####
# -Velmeshev et al.: First tab has the DEGs by cluster

fn_velm <- here(dir_rdata, "Velmeshev_ASD", "aav8130_data-s4_Velmeshev_ASD.xls")

data_velm <- readxl::read_excel(fn_velm, sheet = "ASD_DEGs") |>
  mutate(`Cell type` = factor(`Cell type`))

# `Cell type` `gene ID` Gene …¹ Gene …² Fold …³ Sampl…⁴ q val…⁵ corre…⁶ Epile…⁷ gene …⁸ SFARI…⁹ Satte…˟ Sanders cell …˟
# <chr>       <chr>     <chr>   <chr>     <dbl>   <dbl>   <dbl> <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>  
# 1 L2/3        ENSG0000… TTF2    protei…  -0.226  -0.844 3.87e-5 0.2159… yes     group I no      no      no      no     
# 2 L2/3        ENSG0000… MX2     protei…  -0.159  -1.16  7.00e-5 -0.034… no      group I no      no      no      no     
# 3 L2/3        ENSG0000… ASCC1   protei…  -0.170  -0.428 1.17e-4 0.3560… no      group I no      no      no      no     
# 4 L2/3        ENSG0000… GLRA3   protei…   0.246   0.397 1.53e-4 0.5133… no      group V no      no      no      no 

data_velm |> filter(`q value` < 0.05) |> count(`Cell type`)
# `Cell type`          n
# <chr>            <int>
# 1 AST-FB              15
# 2 AST-PP              60
# 3 Endothelial         60
# 4 IN-PV               28
# 5 IN-SST              29
# 6 IN-SV2C             32
# 7 IN-VIP              82
# 8 L2/3               109
# 9 L4                 108
# 10 L5/6                 8
# 11 L5/6-CC             19
# 12 Microglia           80
# 13 Neu-mat              3
# 14 Neu-NRGN-I          14
# 15 Neu-NRGN-II         26
# 16 Oligodendrocytes    16
# 17 OPC                  3



## extract gene list
genes_velm <- data_velm |> 
  filter(`q value` < 0.05)  |> 
  group_by(`Cell type`) |> 
  group_map(~pull(.x,`gene ID`))

names(genes_velm) <- levels(data_velm$`Cell type`)

#### Batiuk ####
# -Batiuk et al.: Table 1 has DEGs by cluster

fn_batiuk <- here(dir_rdata, "Batiuk_SCZ", "Supplementary_Dataset_Tables_1-4.xlsx")

data_batiuk <- readxl::read_excel(fn_batiuk, sheet = 1)

## fix header ...easier than fixing colnames with skip=2
colnames(data_batiuk) <- data_batiuk[2,]
data_batiuk <- tail(data_batiuk, -2)

dim(data_batiuk)
# [1] 462631     26

## add halves together, convert types
data_batiuk <- rbind(data_batiuk[,1:12], data_batiuk[,15:26]) |> 
  type_convert()

# Cell_type     baseMean log2FoldChange  lfcSE  stat        pvalue      padj     Z    Za Gene       CellFrac SampleFrac
# <chr>            <dbl>          <dbl>  <dbl> <dbl>         <dbl>     <dbl> <dbl> <dbl> <chr>         <dbl>      <dbl>
# 1 ID2_LAMP5_CRH    179.           3.69  0.620   5.96 0.00000000259 0.0000320  5.96  4.16 FP236383.1   0.198           1
# 2 ID2_LAMP5_CRH    393.           2.68  0.534   5.01 0.000000537   0.00332    5.01  2.94 AC060765.2   0.439           1
# 3 ID2_LAMP5_CRH     19.2          2.25  0.468   4.82 0.00000145    0.00598    4.82  2.75 HES4         0.0514          1
# 4 ID2_LAMP5_CRH    241.           0.786 0.169   4.67 0.00000307    0.00947    4.67  2.59 CHD5         0.508           1
# 5 ID2_LAMP5_CRH    175.           0.412 0.0908  4.54 0.00000568    0.0131     4.54  2.48 CDC42BPB     0.421           1

data_batiuk |> 
  filter(padj < 0.05) |>
  count(Cell_type) |> 
  arrange(-n)

# Cell_type               n
# <chr>               <int>
#   1 L4_RORB_SCHLAP1_MET   155
# 2 L4_RORB_SCHLAP1        32
# 3 ID2_LAMP5_CRH          16
# 4 ID2_LAMP5              12
# 5 L2_CUX2_LAMP5_PDGFD    12
# 6 PVALB                  12
# 7 SST_TAC3               12
# 8 ID2_PAX6               10
# 9 VIP                     8
# 10 L5_6_THEMIS_SEMA3A      6
# … with 23 more rows


#### Export Gene lists ####
sn_DEGs <- list(Velmeshev = genes_velm)

save(sn_DEGs, file = here(dir_rdata, "singleNuc_DEGs.Rdata"))



