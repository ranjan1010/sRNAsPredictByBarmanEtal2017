# this functions helps us remove previously created workspace variables....
rm(list=ls());

# this library used to handle string related functions.......
library("stringr");

# use to get input complete genome sequence data in .fasta file from user 
inputFastaFile <- readline("please enter text complete genome sequence data file in .fasta file: ");

#Read text data with all lines 
readFastaFile <- paste0(readLines(inputFastaFile));
readFastaFile

# use to get input start end location information for a complete sequence in .csv file from user 
inputIndexFile <- readline("please enter source file name which consists start end location of sequence in extension (.csv): ");
indexFile <- read.csv(inputIndexFile, header = FALSE);
indexFile

length(readFastaFile)

# This function used to concatenated all lines string into one line. This function will take only read parameter(read text file name)  
concatAllLinesInOne <- function(readFile)
{
  finalString = " ";
  # the first line for any fasta file is fasta header file that's why loop started from line 2 
  for(i in 2 : length(readFile))
  {
    finalString = paste0(finalString,readFile[i])
  }
  trimFinalString <- str_trim(finalString);
  return(trimFinalString);
}

singleLineString = concatAllLinesInOne(readFastaFile);
singleLineString

# This is use to temporary  store the data

tempMatrix = matrix(0, nrow = length(indexFile[,1]), ncol = 5);

# Header initialization
tempMatrix[1,1] = toString(indexFile[1,1]);
tempMatrix[1,2] = toString(indexFile[1,2]);
tempMatrix[1,3] = toString(indexFile[1,3]);
tempMatrix[1,4] = toString(indexFile[1,4]);
tempMatrix[1,5] = "Sequence";
indexFile[1,1]

# Iterate for all index file........ 
for (i in 2 : length(indexFile[,1]))
{
  
  # All previous data were stored also in new matrix.
  tempMatrix[i,1] = toString(indexFile[i,1]);
  tempMatrix[i,2] = toString(indexFile[i,2]);
  tempMatrix[i,3] = toString(indexFile[i,3]);
  tempMatrix[i,4] = toString(indexFile[i,4]);
  
  # This section will store sRNA sequence
  tempMatrix[i,5] = toString(substr(singleLineString, as.numeric(as.numeric_version(indexFile[i,2])), as.numeric(as.numeric_version(indexFile[i,3]))));
  
  
}
tempMatrix

outputFile <- readline("please enter output file name for sRNA data in .csv file: ");
write.table(tempMatrix,file = outputFile, sep=",",quote = FALSE, col.names=FALSE, row.names = FALSE);


