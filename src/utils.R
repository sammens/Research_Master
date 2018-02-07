#################################################################
#  Contributor: Samuel Ofosu Mensah
#################################################################

### A function to plot multiple histograms
histplot = function(x, r, co, st){
  # where x is the input is a matrix and is the new data after we have removed the unwanted features
  # r is the number of graphs we want on a row
  # co is the number of graphs we want on a column
  # st is which index to start
  # sp is which index to stop
  x = t(x) 	#transpose data
  par(mfrow = c(r,co))
  for(i in st:(st + ((r*co)-1))){
    graphs = hist(x[,i], xlab = paste(colnames(x)[i]), main = paste('Histogram of ', colnames(x)[i]))
  }
  graphs
}

