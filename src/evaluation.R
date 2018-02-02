###############################################################
#  Contributor: Samuel Ofosu Mensah
###############################################################

## Applying the SAMseq function

library(samr)
samfit = SAMseq(x = data$x, y = data$y, resp.type = "Two class unpaired", geneid = data$geneid, genenames = data$genenames)
