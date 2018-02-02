###############################################################
#  Contributor: Samuel Ofosu Mensah
###############################################################

## Applying the SAMseq function

library(samr)
trace('samr.estimate.depth', edit = TRUE)
samfit = SAMseq(x = mydata$x, y = mydata$y, resp.type = "Two class unpaired", geneid = mydata$geneid, genenames = mydata$genenames)
