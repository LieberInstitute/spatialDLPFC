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
k = c(3,9,16)
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
  if(k[i]==9){
    mycolors <- mycolors[c(1,2,3,5,7,6,9,4,8)]
  }
  if(k[i]==16){
    mycolors <- mycolors[c(13,11,6,15,4,16,7,12,2,1,14,5,9,8,3,10)]
  }
  png(paste0("/Users/abbyspangler/Desktop/colorbar",k[i],".png"))
  color.bar(mycolors = mycolors)
  dev.off()
}

