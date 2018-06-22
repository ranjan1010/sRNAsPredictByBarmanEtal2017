# this functions helps us remove previously created workspace variables....
rm(list=ls());

# this library used to handle string related functions.......
library("stringr");

# use to get input start and end position information of all CDS for in .csv file format from user 
inputCodingIndexFile <- readline("please enter source file name which consists start and end position of all CDS in extension (.csv): ");
CodingIndexFile <- read.csv(inputCodingIndexFile, header = TRUE);
CodingIndexFile

# use to get end position of complete genome of particular strain 
getLastPositionOfCompleteGemone <- readline("please enter the end position of complete genome of particular strain : ");
lastPositionOfCompleteGemone <- as.numeric(getLastPositionOfCompleteGemone);

# find the column index of start position of CDS
columnIndexOfStart <- match(toString("Start"),names(CodingIndexFile));
columnIndexOfStart

# find the column index of stop position of CDS
columnIndexOfStop <- match(toString("Stop"),names(CodingIndexFile));
columnIndexOfStop


tempMatrix = matrix(0, nrow = length(CodingIndexFile[,1]) + 1, ncol = 4);

# Header initialization
tempMatrix[1,1] = "IntergenicRegion";
tempMatrix[1,2] = "Start";
tempMatrix[1,3] = "Stop";
tempMatrix[1,4] = "Length";

tempCountForNonIntergenic = 0;
tempCountForIntergenic = 1;
# Iterate for all index file........ 
for (i in 1 : length(CodingIndexFile[,1]))
{
  if(i == 1)
  {
    tempMatrix[i+1,1] = paste0("InterGenicRegion",toString(as.numeric(i)));
    tempMatrix[i+1,2] = 1;
    tempMatrix[i+1,3] = as.numeric(CodingIndexFile[i,columnIndexOfStart]) - 1;
    tempMatrix[i+1,4] = as.numeric(CodingIndexFile[i,columnIndexOfStart]) - 1;
  }
  else if(i == length(CodingIndexFile[,1]))
  {
    tempMatrix[i+1,1] = paste0("InterGenicRegion",toString(as.numeric(tempCountForIntergenic)+1));
    tempMatrix[i+1,2] = as.numeric(CodingIndexFile[i,columnIndexOfStop]) + 1;
    tempMatrix[i+1,3] = lastPositionOfCompleteGemone;
    tempMatrix[i+1,4] = lastPositionOfCompleteGemone - as.numeric(CodingIndexFile[i,columnIndexOfStop]);
  }
  else
  {
    if(as.numeric(CodingIndexFile[i-1,columnIndexOfStop]) < as.numeric(CodingIndexFile[i,columnIndexOfStart]))
    {
      tempCountForIntergenic = tempCountForIntergenic + 1
      tempMatrix[i+1,1] = paste0("InterGenicRegion",toString(as.numeric(tempCountForIntergenic)));
      tempMatrix[i+1,2] = as.numeric(CodingIndexFile[i-1,columnIndexOfStop]) + 1;
      tempMatrix[i+1,3] = as.numeric(CodingIndexFile[i,columnIndexOfStart]) - 1;
      tempMatrix[i+1,4] = (as.numeric(CodingIndexFile[i,columnIndexOfStart]) - as.numeric(CodingIndexFile[i-1,columnIndexOfStop])) - 1;
    }
    else
    {
      tempCountForNonIntergenic = tempCountForNonIntergenic + 1;
      tempMatrix[i+1,1] = paste0("Non_InterGenicRegion",toString(as.numeric(tempCountForNonIntergenic)));
    }
    
  }
}
tempMatrix

outputFile <- readline("please enter output file name for intergenic regions index in .csv file: ");
write.table(tempMatrix,file = outputFile, sep=",",quote = FALSE, col.names=FALSE, row.names = FALSE);


