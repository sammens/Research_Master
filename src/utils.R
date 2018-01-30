"""Contributor: Samuel Ofosu Mensah"""

"""A function to plot multiple histograms"""
histplot = function(x, r, co, st){
     #where x is the input is a matrix and is the new data after we have removed the unwanted features
     #r is the number of graphs we want on a row
     #co is the number of graphs we want on a column
     #st is which index to start
     #sp is which index to stop
     par(mfrow = c(r,co))
     for(i in st:(st + ((r*co)-1)){
         graphs = hist(x[,i])
     }
     graphs
 
}

