library(RColorBrewer)


## Make the plot
mycolors <- brewer.pal(8, "Dark2")

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
    ylab = ""#,
    #breaks = breaks,
    #nlevel = length(mypal),
    #axis.args = axis.args
  )
  return(bar)
}

png("/Users/abbyspangler/Desktop/test_bar.png")
color.bar(mycolors = mycolors)
dev.off()
