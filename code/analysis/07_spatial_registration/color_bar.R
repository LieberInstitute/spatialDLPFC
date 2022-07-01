library(RColorBrewer)

color.bar <- function(mycolors){
  n.color <- length(mycolors)
  #par(mar = mar) #adjust width and height with par 
  bar<-image(
    x = 1:n.color,
    y = 1,
    z = as.matrix(1:n.color),
    col = mycolors,
    xaxt = "n",
    yaxt = "n",
    xlab = "",
    ylab = ""
  )
  return(bar)
}

## Make the plot
k = c(3,7,9,16,28)
for(i in seq_along(k) ){
  print(k[i])
  mycolors <- Polychrome::palette36.colors(k[i])
  names(mycolors)<- c(1:k[i])
  print(mycolors)
  #print plot
  # png(paste0("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/plots/07_spatial_registration_",k[i],".png"))
  # color.bar(mycolors = mycolors)
  # dev.off()
  # 
  if(k[i]==7){
    mycolors <- mycolors[c(1,5,7,2,3,4,6)]
  }
  if(k[i]==9){
    mycolors <- mycolors[c(1,2,6,9,4,7,3,5,8)]
  }
  if(k[i]==16){
    mycolors <- mycolors[c(1,2,14,6,15,11,13,4,16,7,12,5,9,3,8,10)]
  }
  if(k[i]==28){
    mycolors <- mycolors[c(16,17,20,5,28,4,13,19,7,14,11,26,27,23,12,15,22,25,5,10,24,8,9,1,2,3)]
  }
  png(paste0("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/plots/07_spatial_registration/colobar_k",k[i],".png"))
  color.bar(mycolors = mycolors)
  dev.off()
}

