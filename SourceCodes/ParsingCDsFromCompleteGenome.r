# this functions helps us remove previously created workspace variables....
rm(list=ls());

# this library used to handle string related functions.......
library("stringr");

# use to get complete list of all CDs file in .txt file from user 
inputTextfile <- readline("please enter the complete list of all CDs file in .txt file format: ");

#Read text data with all lines 
readTextFile <- paste0(readLines(inputTextfile));
readTextFile

# This section will find all the line numbers with keyword "CDS" and return the particular line for CDs
allCDs <- grep("CDS",readTextFile,value = TRUE);
allCDs

# This section will find all the forward CDS
ForwardCDs <- grep("complement",allCDs,value = TRUE,invert = TRUE);
ForwardCDs

# This section will remove the CDS and space from string
ForwardCDs1 <- gsub("     CDS             ","",ForwardCDs);

# This section will replace .. to , 
ForwardCDs2 <- gsub("\\.\\.",",",ForwardCDs1);
ForwardCDs2

outputFile <- readline("please enter the output file name for CDS start and end position of coding part in .csv file format: ");
write.table(ForwardCDs2,file = outputFile, sep=",",quote = FALSE, col.names=FALSE, row.names = FALSE);
