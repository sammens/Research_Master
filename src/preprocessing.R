##############################################################
#  Contributor: Samuel Ofosu Mensah
##############################################################

### A function to clean TCGA dataset
clean_data = function(x){
  # input is unprocessed data from TCGA
  BrcaHeader = names(x)		#names of data
  Patient_Id = BrcaHeader[1]	#Patient ID
  BrcaHeader = BrcaHeader[-1]	#gene names
  BrcaHeader = unlist(strsplit(BrcaHeader, '|', fixed=T))  #splitting characters on either sides of '|' in the header
  BrcaHeader = BrcaHeader[seq(1,length(BrcaHeader),2)]	#removing numbers
  names(x) = c(Patient_Id, BrcaHeader)	#changing the header names of the data
  x = x[, !grepl("\\?", names(x))]	#removing all genes starting with "?"
  return(x)
}

### A function to label samples as tumour or normal
sampletwt = function(x){
  # input is processed data from TCGA
  s <- as.numeric(substr(unlist(strsplit(x, '-', fixed = TRUE))[4],1,2)) #function method
  #if statement to select tumour, normal, and control
  # tumour: 1 to 9 
  # normal: 10 to 19
  # control: 20 to 29
  if(s %in% 1:9) {
    return("tumour")
  } else if(s %in% 10:19) {
      return("normal")
  } else {
      return("control")
    }
}

### A function to apply "sampletwt"
label.data = function(x){
  # input is processed data from TCGA  
  catpatients = apply(x[1], 1, sampletwt)	#list of all categories
  x = as.data.frame(append(x, list(groups=catpatients), after = 1))   #adding 'catpatients' to x
  return(x)
}


### A function for Coefficient of Variation
CV = function(x){
  # where x is the input is a matrix
  # CV stands for Coefficient of Variation
  cv = (sd(x)/mean(x))
  return(cv)
}

### A function to remove all columns with mean == 0
remove_NaN = function(x){
  #where x is the input is a matrix
  x = data.frame(x)
  cv = apply(x,2,CV)	   	#applying the CV(x) function
  x = x[,-which(cv == 'NaN')] 	#remove all features with mean == 0
  return(x)
}

### A function for Miller's test statistic
Miller = function(x){
  #where x is the input is a matrix
  c_0 = 1/3	 		#c_0 to represent population Coefficient of Variation
  m = length(x) - 1		#length of sample
  cv = CV(x) 			#using the Coefficient of Variation function
  z = ((sqrt(m)) * (cv - c_0))/(sqrt(c_0 * (0.5 + c_0**2))) #Miller's statistic
  decision = (1 - pnorm(z)) > 0.05	#deciding null hypothesis 
  return(decision)
}


### A function to round-up real numbers to integers
rnd = function(x){
  x = trunc(x + sign(x) * 0.5)
  return(x)
}

### A function to generate 'x' data for SAMseq
gen_xdata = function(x){
  x = remove_NaN(x)
  Miller.test = apply(x,2,Miller)
  x = x[,which(Miller.test == TRUE)]
  x = apply(x,2,rnd)
  mode(x) = 'integer'
  x = t(x)
  return(x)
}

### A function to generate 'y' data for SAMseq
gen_ydata = function(x){
  x = x$groups
  x = as.character(x)
  x[x == 'tumour'] = 1
  x[x == 'normal'] = 2
  x = as.numeric(x)
  return(x)
}

### A function to combine the required data
list.data = function(x){
  xx = x[,-c(1,2)]
  xx = gen_xdata(xx)
  geneid = as.character(1:nrow(xx))
  genenames = row.names(xx)
  data = list(x = xx, y = gen_ydata(x), geneid = geneid, genenames = genenames)
  return(data)
}

### A function to return all written functions
all.func = function(x){
  #using remove_NaN to remove all columns with mean == 0
  x = remove_NaN(x)
  #applying Miller's test after removing coloumns with mean == 0
  mil = apply(x,2,Miller)
  return(list('data' = x, 'Miller' = mil))
}
