"""Contributor: Samuel Ofosu Mensah"""

"""A function for Coefficient of Variation"""
CV = function(x){
  #input is a matrix
  #sd for standard deviation
  cv = (sd(x)/mean(x))
  return(cv)
}

"""A function to remove all columns with mean == 0"""
remove_NaN = function(x){
    #input is a matrix
    cv = apply(x,2,CV)	#applying the CV(x) function
    x = x[,-which(cv == 'NaN')] 	#remove all features with mean == 0
    return(x)
}

"""A function for Miller's test statistic"""
Miller = function(x){
    #input is a matrix
    c_0 = 1/3	 	#c_0 to represent population Coefficient of Variation
    m = length(x) - 1	#length of sample
    cv = CV(x) 		#using the Coefficient of Variation function
    z = ((sqrt(m)) * (cv - c_0))/(sqrt(c_0 * (0.5 + c_0**2))) #Miller's statistic
    decision = z > 1.64	#deciding null hypothesis 
    return(decision)
}
