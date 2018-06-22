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
  inputIntergenicRegionsFile <- paste0(toString(intergenicRegionsSequenceFile[IGR,1]),".csv");
  intergenicRegionsFile <- read.csv(inputIntergenicRegionsFile, header = TRUE);
  
  # Find the column index of intergenic regions sequence
  matchedColumnName = match("Sequence",names(intergenicRegionsFile));
  matchedColumnName
  
  # This function used to calculate Oligonucleotide frequency of a prticular sequence
  calculateOligonucleotideFrequency <- function(sequence)
  {
    
    # Initialized the all oligonucleotide. The result will be come out in "A","T","G","C" order 
    oligonucleotide <- c("A","T","G","C");
    
    # This matrix use to store Oligonucleotide frequency with Oligonucleotide
    tempMatrixForOligonucleotide <- matrix(0, nrow = 1, ncol = length(oligonucleotide));
    for(i in 1 : length(oligonucleotide))
    {
      # This is use to calculate oligonucleotide frequency in "A","T","G","C" order 
      tempMatrixForOligonucleotide[1,i] = str_count(sequence, oligonucleotide[i])/str_length(sequence);
    }
    tempMatrixForOligonucleotide
    return(tempMatrixForOligonucleotide);
  }
  
  # This function used to calculate dinucleotide frequency of a particular sequence
  calculateDinucleotideFrequency <- function(sequence)
  {
    
    # Initialized the all dinucleotide The result will be come out in "AA","AT","AG","AC","TA","TT","TG","TC","GA","GT","GG","GC","CA","CT","CG","CC" order 
    dinucleotide <- c("AA","AT","AG","AC","TA","TT","TG","TC","GA","GT","GG","GC","CA","CT","CG","CC");
    
    # This matrix use to store dinucleotide frequency with dinucleotide
    tempMatrixForDinucleotide <- matrix(0, nrow = 1, ncol = length(dinucleotide));
    for(i in 1 : length(dinucleotide))
    {
      # This Count variable use to count the particular dinucleotide
      Count = 0;
      for(j in 1 : str_length(sequence))
      {
        # This is use to count all dinucleotide with window size 2 for this reason all possible window lenght 2 (j, j+1) will consider here  
        tempCount = str_count(substring(sequence, j,j+1), dinucleotide[i]);
        Count = Count + tempCount;
      }
      Count
      # This is use to calculate Digonucleotide frequency in "AA","AT","AG","AC","TA","TT","TG","TC","GA","GT","GG","GC","CA","CT","CG","CC" order 
      tempMatrixForDinucleotide[1,i] = Count/(str_length(sequence) - 1) ;
    }
    tempMatrixForDinucleotide
    return(tempMatrixForDinucleotide);
  }
  
  # This function used to calculate trinucleotide frequency of a particular sequence
  calculateTrinucleotideFrequency <- function(sequence)
  {
    
    # Initialized the all trinucleotide The result will be come out the following oder 
    trinucleotide <- c("AAA","ATA","AGA","ACA","AAT","ATT","AGT","ACT","AAG","ATG","AGG","ACG",
                       "AAC","ATC","AGC","ACC","TAA","TTA","TGA","TCA","TAT","TTT","TGT","TCT",
                       "TAG","TTG","TGG","TCG","TAC","TTC","TGC","TCC","GAA","GTA","GGA","GCA",
                       "GAT","GTT","GGT","GCT","GAG","GTG","GGG","GCG","GAC","GTC","GGC","GCC",
                       "CAA","CTA","CGA","CCA","CAT","CTT","CGT","CCT","CAG","CTG","CGG","CCG",
                       "CAC","CTC","CGC","CCC");
    
    # This matrix use to store trigonucleotide frequency with trinucleotide
    tempMatrixForTrinucleotide <- matrix(0, nrow = 1, ncol = length(trinucleotide));
    for(i in 1 : length(trinucleotide))
    {
      # This Count variable use to count the particular trinucleotide
      Count = 0;
      for(j in 1 : str_length(sequence))
      {
        # This is use to count all trinucleotide with window size 3 for this reason all possible window lenght 3 (j, j+2) will consider here  
        tempCount = str_count(substring(sequence, j,j+2), trinucleotide[i]);
        Count = Count + tempCount;
      }
      Count
      # This is use to calculate trigonucleotide frequency in "AAA","ATA","AG"A,"ACA" etc. ......... order 
      tempMatrixForTrinucleotide[1,i] = Count/(str_length(sequence) - 2) ;
    }
    tempMatrixForTrinucleotide
    return(tempMatrixForTrinucleotide);
  }
  
  tempMatrix <- matrix(0, nrow = length(intergenicRegionsFile[,1]) + 1, ncol = 88);
  
  # this section use to header initialization for output matrix 
  for(i in 1 : length(names(intergenicRegionsFile)) - 1)
  {
    tempMatrix[1,i] = names(intergenicRegionsFile[i]);
  }
  
  headerForOligonucleotide <- c("A","T","G","C")
  
  # this section will initialize the Oligonucleotide where 4 (input file header name) 
  # + 4 (oligonuclotide frequency)
  for(i in length(names(intergenicRegionsFile)) : 8)
  {
    tempMatrix [1,i] = headerForOligonucleotide[i -(length(names(intergenicRegionsFile)) - 1)] ;
  }
  
  
  headerForDinucleotide <- c("AA","AT","AG","AC","TA","TT","TG","TC","GA","GT","GG","GC","CA","CT","CG","CC");
  
  # this section will initialize the Oligonucleotide where 4 (input file header name)  
  #   + 4 (oligonuclotide frequency) + 16 (dinucleotide frequency)
  for(i in 9 : 24)
  {
    tempMatrix [1,i] = headerForDinucleotide[i - 8];
  }
  
  headerForTrinucleotide <- c("AAA","ATA","AGA","ACA","AAT","ATT","AGT","ACT","AAG","ATG","AGG","ACG",
                              "AAC","ATC","AGC","ACC","TAA","TTA","TGA","TCA","TAT","TTT","TGT","TCT",
                              "TAG","TTG","TGG","TCG","TAC","TTC","TGC","TCC","GAA","GTA","GGA","GCA",
                              "GAT","GTT","GGT","GCT","GAG","GTG","GGG","GCG","GAC","GTC","GGC","GCC",
                              "CAA","CTA","CGA","CCA","CAT","CTT","CGT","CCT","CAG","CTG","CGG","CCG",
                              "CAC","CTC","CGC","CCC");
  
  # this section will initialize the trinucleotide frquency where 4 (input file header name) 
  # + 4 (oligonuclotide frequency) + 16 (dinucleotide frequency) + 64 (trinucleotide fequency)
  for(i in 25 : 88)
  {
    tempMatrix [1,i] = headerForTrinucleotide[i - 24];
  }
  
  # Iterate for all rows for input file
  for(i in 1 : length(intergenicRegionsFile[,1]))
  {
    tempMatrix[i+1,1] = toString(intergenicRegionsFile[i,1]);
    tempMatrix[i+1,2] = toString(intergenicRegionsFile[i,2]);
    tempMatrix[i+1,3] = toString(intergenicRegionsFile[i,3]);
    tempMatrix[i+1,4] = toString(intergenicRegionsFile[i,4]);
    #----------------------------------------------- oligonucleotide frequency calculation start here ----------------------
    #This function will calculate oligonucleotide frequency for a given sequence
    tempOligonucleotideFrequency = calculateOligonucleotideFrequency(intergenicRegionsFile[i,matchedColumnName]);
    
    # This section use to store the oligonucleotide frequency value in output matrix
    for(j in 1: length(tempOligonucleotideFrequency[1,]))
    {
      tempMatrix[i+1,j+4] = tempOligonucleotideFrequency[1,j];
    }
    
    #----------------------------------------------- oligonucleotide frequency calculation end here ----------------------
    
    #----------------------------------------------- dinucleotide frequency calculation start here ----------------------
    
    #This function will calculate dinucleotide frequency for a given sequence
    tempDinucleotideFrequency = calculateDinucleotideFrequency(intergenicRegionsFile[i,matchedColumnName]);
    
    # This section use to store the dinucleotide frequency value in output matrix
    for(j in 1: length(tempDinucleotideFrequency[1,]))
    {
      tempMatrix[i+1,j+8] = tempDinucleotideFrequency[1,j];
    }
    
    #----------------------------------------------- dinucleotide frequency calculation end here ----------------------
    
    #----------------------------------------------- trinucleotide frequency calculation start here ----------------------
    
    #This function will calculate trinucleotide frequency for a given sequence
    tempTrinucleotideFrequency = calculateTrinucleotideFrequency(intergenicRegionsFile[i,matchedColumnName]);
    
    # This section use to store the trinucleotide frequency value in output matrix
    for(j in 1: length(tempTrinucleotideFrequency[1,]))
    {
      tempMatrix[i+1,j+24] = tempTrinucleotideFrequency[1,j];
    }
    
    #----------------------------------------------- trinucleotide frequency calculation end here ----------------------
  }
  tempMatrix
  outputFile <- paste0("AllNucleotideFrequencyOf",toString(intergenicRegionsSequenceFile[IGR,1]),".csv");
  write.table(tempMatrix,file = outputFile, sep=",",quote = FALSE, col.names=FALSE, row.names = FALSE);
  
}
