# this functions helps us remove previously created workspace variables....
rm(list=ls());

# this library used to handle string related functions.......
library("stringr");

# use to get input intergenic regions sequence in .csv file from user 
inputIntergenicRegionsSequenceFile <- readline("please enter intergenic regions sequence file in extension (.csv): ");
intergenicRegionsSequenceFile <- read.csv(inputIntergenicRegionsSequenceFile, header = FALSE);
intergenicRegionsSequenceFile

# Use to get the length of nucleotide for sliding window  
inputSildingWindowSize <- readline("please enter the sliding window size: ");
sildingWindowSize <- as.numeric(inputSildingWindowSize);
sildingWindowSize

# Use to get the step size for sliding window
inputStepSize <- readline("please enter the step size for sliding window: ");
stepSize <- as.numeric(inputStepSize);
stepSize


for(i in 2 : length(intergenicRegionsSequenceFile[,1]))
{
  if((as.numeric(as.numeric_version(as.numeric(as.numeric_version(intergenicRegionsSequenceFile[i,4])))) >= 50) && (as.numeric(as.numeric_version(as.numeric(as.numeric_version(intergenicRegionsSequenceFile[i,4])))) <= sildingWindowSize))
  {
    tempMatrix <- matrix(0,nrow = 2,ncol = 5);
    
    tempMatrix[1,1] = "IntergenicRegionName";
    tempMatrix[1,2] = "WindowName";
    tempMatrix[1,3] = "WindowStart";
    tempMatrix[1,4] = "WindowStop";
    tempMatrix[1,5] = "Sequence";
    
    tempMatrix[2,1] = toString(intergenicRegionsSequenceFile[i,1]);
    tempMatrix[2,2] = "Window1";
    tempMatrix[2,3] = 1
    tempMatrix[2,4] = as.numeric(as.numeric_version(intergenicRegionsSequenceFile[i,4]));
    tempMatrix[2,5] = toString(substr(intergenicRegionsSequenceFile[i,5],1,as.numeric(as.numeric_version(intergenicRegionsSequenceFile[i,4]))));
    tempMatrix
    
    outputFile <- paste0(toString(intergenicRegionsSequenceFile[i,1]),".csv");
    write.table(tempMatrix,file = outputFile, sep=",",quote = FALSE, col.names=FALSE, row.names = FALSE);
  }
  else
  {
    noOfWindow = 0;
    noOfWindow <- ceiling((as.numeric(as.numeric_version(intergenicRegionsSequenceFile[i,4])) - sildingWindowSize)/stepSize) + 1;
    noOfWindow
    if(noOfWindow > 1)
    {
      tempMatrix <- matrix(0, nrow = noOfWindow + 1, ncol = 5);
      tempMatrix[1,1] = "IntergenicRegionName";
      tempMatrix[1,2] = "WindowName";
      tempMatrix[1,3] = "WindowStart";
      tempMatrix[1,4] = "WindowStop";
      tempMatrix[1,5] = "Sequence";
        
      for(j in 1 : noOfWindow)
      {
        tempMatrix[j+1,1] = toString(intergenicRegionsSequenceFile[i,1]);
        tempMatrix[j+1,2] = paste0("Window",toString(j));
        if(j == 1)
        {
          tempMatrix[j+1,3] = 1;
          tempMatrix[j+1,4] = sildingWindowSize;
          tempMatrix[j+1,5] = toString(substr(intergenicRegionsSequenceFile[i,5],1,sildingWindowSize));
        }
        else if(j == noOfWindow)
        {
          tempMatrix[j+1,3] = (stepSize * (j -1)) + 1;
          tempMatrix[j+1,4] = as.numeric(as.numeric_version(intergenicRegionsSequenceFile[i,4]));
          tempMatrix[j+1,5] = toString(substr(intergenicRegionsSequenceFile[i,5],((stepSize * (j -1)) + 1),as.numeric(as.numeric_version(intergenicRegionsSequenceFile[i,4]))));
        }
        else
        {
          tempMatrix[j+1,3] = (stepSize * (j -1)) + 1;
          tempMatrix[j+1,4] = (stepSize * (j -1)) + sildingWindowSize ;
          tempMatrix[j+1,5] = toString(substr(intergenicRegionsSequenceFile[i,5],((stepSize * (j -1)) + 1), ((stepSize * (j -1)) + sildingWindowSize)));
        }
      }
      tempMatrix
      
      outputFile <- paste0(toString(intergenicRegionsSequenceFile[i,1]),".csv");
      write.table(tempMatrix,file = outputFile, sep=",",quote = FALSE, col.names=FALSE, row.names = FALSE);
    }
  }
}