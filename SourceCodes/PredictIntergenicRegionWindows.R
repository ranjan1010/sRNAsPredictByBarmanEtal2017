# this functions helps us remove previously created workspace variables....
rm(list=ls());

# this library used to handle string related functions.......
library("stringr");

# use to get input intergenic regions sequence in .csv file from user 
inputIntergenicRegionsSequenceFile <- readline("please enter intergenic regions sequence file in extension (.csv): ");
intergenicRegionsSequenceFile <- read.csv(inputIntergenicRegionsSequenceFile, header = FALSE);
intergenicRegionsSequenceFile

for(IGR in 2 : length(intergenicRegionsSequenceFile[,1]))
{
  # This section will classify the data
  system(paste0("svm_classify.exe"," ",toString(intergenicRegionsSequenceFile[IGR,1]),"TriNucleotideCompositionSVMData",".txt"," ",
                "ModelTriNucleotideCompositionOneistoTwo7_7_4_3", " ","PredictionScoreOf",toString(intergenicRegionsSequenceFile[IGR,1]),".txt"));
}