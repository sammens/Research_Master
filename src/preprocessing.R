"""Contributor: Samuel Ofosu Mensah"""

"""A function for Coefficient of Variation"""
CV = function(x){
  #where x is the input is a matrix
  #sd for standard deviation
  cv = (sd(x)/mean(x))
  return(cv)
}

"""A function to remove all columns with mean == 0"""
remove_NaN = function(x){
    #where x is the input is a matrix
    cv = apply(x,2,CV)	#applying the CV(x) function
    x = x[,-which(cv == 'NaN')] 	#remove all features with mean == 0
    return(x)
}

"""A function for Miller's test statistic"""
Miller = function(x){
    #where x is the input is a matrix
    c_0 = 1/3	 	#c_0 to represent population Coefficient of Variation
    m = length(x) - 1	#length of sample
    cv = CV(x) 		#using the Coefficient of Variation function
    z = ((sqrt(m)) * (cv - c_0))/(sqrt(c_0 * (0.5 + c_0**2))) #Miller's statistic
    decision = (1 - pnorm(z)) > 0.05	#deciding null hypothesis 
    return(decision)
}

"""A function to plot multiple histograms for the selected genes"""
histplot = function(x, r, co, st){
    #where x is the input is a matrix and is the new data after we have removed the unwanted features
    #r is the number of graphs we want on a row
    #co is the number of graphs we want on a column
    #st is which index to start
    #sp is which index to stop
    par(mfrow = c(r,co))
    for(i in st:(st + 19)){
        graphs = hist(x[,i])
    }
    graphs
}

"""A function to return all written functions"""
all.func = function(x){
    #using remove_NaN to remove all columns with mean == 0
    x = remove_NaN(x)
    #applying Miller's test after removing coloumns with mean == 0
    mil = apply(x,2,Miller)
    return(list('data' = x, 'Miller' = mil))
}
