# this functions helps us remove previously created workspace variables....
rm(list=ls());

# this library used to handle string related functions.......
library("stringr");

# use to get input intergenic regions sequence in .csv file from user 
inputIntergenicRegionsSequenceFile <- readline("please enter intergenic regions sequence file in extension (.csv): ");
intergenicRegionsSequenceFile <- read.csv(inputIntergenicRegionsSequenceFile, header = FALSE);
intergenicRegionsSequenceFile

# This boolean variable decide whether user want to generate mono-nucleotide features. if yes supply 1 otherwise supply 0 for no
inputMonoNucleotideDicision <- readline("please enter 1 if user want to generate mono-nucleotide features, otherwise enter 0 : ");
isMonoNucleotide <- as.numeric(inputMonoNucleotideDicision);
isMonoNucleotide

# This boolean variable decide whether user want to generate di-nucleotide features. if yes supply 1 otherwise supply 0 for no
inputDiNucleotideDicision <- readline("please enter 1 if user want to generate di-nucleotide features, otherwise enter 0 : ");
isDiNucleotide <- as.numeric(inputDiNucleotideDicision);
isDiNucleotide

# This boolean variable decide whether user want to generate tri-nucleotide features. if yes supply 1 otherwise supply 0 for no
inputTriNucleotideDicision <- readline("please enter 1 if user want to generate tri-nucleotide features, otherwise enter 0 : ");
isTriNucleotide <- as.numeric(inputTriNucleotideDicision);
isTriNucleotide

for(IGR in 2 : length(intergenicRegionsSequenceFile[,1]))
{
  
  # use to get input features of intergenic region in .csv from user 
  inputsIntergenicRegionFeaturesFile <- paste0("AllNucleotideFrequencyOf",toString(intergenicRegionsSequenceFile[IGR,1]),".csv");
  intergenicRegionFeaturesFile <- read.csv(inputsIntergenicRegionFeaturesFile, header = TRUE);
  intergenicRegionFeaturesFile
	
  
  # this is use to get output file name from user
	outputFileNameFromUser <- toString(intergenicRegionsSequenceFile[IGR,1]);

	# Initialized the all mono-nucleotide.
	mononucleotide <- c("A","T","G","C");

	# Initialized the all di-nucleotide.
	dinucleotide <- c("AA","AT","AG","AC","TA","TT","TG","TC","GA","GT","GG","GC","CA","CT","CG","CC");

	# Initialized the all tri-nucleotide.
	trinucleotide <- c("AAA","ATA","AGA","ACA","AAT","ATT","AGT","ACT","AAG","ATG","AGG","ACG",
					   "AAC","ATC","AGC","ACC","TAA","TTA","TGA","TCA","TAT","TTT","TGT","TCT",
					   "TAG","TTG","TGG","TCG","TAC","TTC","TGC","TCC","GAA","GTA","GGA","GCA",
					   "GAT","GTT","GGT","GCT","GAG","GTG","GGG","GCG","GAC","GTC","GGC","GCC",
					   "CAA","CTA","CGA","CCA","CAT","CTT","CGT","CCT","CAG","CTG","CGG","CCG",
					   "CAC","CTC","CGC","CCC");
	  
	  if( isMonoNucleotide && isDiNucleotide && isTriNucleotide)
	  {
		# This section will generate all nucleotide composition (84) features along with indicator of class (+1 or -1)
		tempMatrixAllNucleotideComposition = matrix(0, nrow = length(intergenicRegionFeaturesFile[,1]), ncol = length(mononucleotide) + 
											 length(dinucleotide) + length(trinucleotide) + 1);
		
		# --------------------------------- Features vector generate for intergenicRegion set start here ------------------------ 
		
		for(i in 1 : length(intergenicRegionFeaturesFile[,1]))
		{
		  tempMatrixAllNucleotideComposition[i,1] = toString("0");
		  
		  # Use to count temp cloumn number
		  tempColumnCount = 2;
		  
		  # -------------------------------------- mono-nucleotide assign here for intergenicRegion set -------------------------------------
		  for(j in 1 : length(mononucleotide))
		  {
			columnNoMono = match(mononucleotide[j],names(intergenicRegionFeaturesFile));
			tempMatrixAllNucleotideComposition[i,tempColumnCount] = paste0(toString(tempColumnCount - 1),":", toString(intergenicRegionFeaturesFile[i,columnNoMono]));
			tempColumnCount = tempColumnCount + 1;
		  }
		  
		  # -------------------------------------- di-nucleotide assign here for intergenicRegion set---------------------------------------
		  for(j in 1 : length(dinucleotide))
		  {
			columnNoDi = match(dinucleotide[j],names(intergenicRegionFeaturesFile));
			tempMatrixAllNucleotideComposition[i,tempColumnCount] = paste0(toString(tempColumnCount - 1),":", toString(intergenicRegionFeaturesFile[i,columnNoDi]));
			tempColumnCount = tempColumnCount + 1;
		  }
		  
		  # -------------------------------------- tri-nucleotide assign here for intergenicRegion set ---------------------------------------
		  for(j in 1 : length(trinucleotide))
		  {
			columnNoTri = match(trinucleotide[j],names(intergenicRegionFeaturesFile));
			tempMatrixAllNucleotideComposition[i,tempColumnCount] = paste0(toString(tempColumnCount - 1),":", toString(intergenicRegionFeaturesFile[i,columnNoTri]));
			tempColumnCount = tempColumnCount + 1;
		  } 
		}
		
		# --------------------------------- Features vector generate for intergenicRegion set end here !!!------------------------ 
		
		
		tempMatrixAllNucleotideComposition
		outputFile <- paste0(toString(outputFileNameFromUser),"AllNucleotideCompositionSVMData", ".txt");
		write.table(as.data.frame(tempMatrixAllNucleotideComposition),file = outputFile,sep=" ",quote = FALSE,col.names=FALSE, row.names = FALSE);
	  
		}else if(isMonoNucleotide && isDiNucleotide){
		  
		  # This section will generate mono and di nucleotide composition (20) features along with indicator of class (+1 or -1)
		  tempMatrixMonoDiNucleotideComposition = matrix(0, nrow = length(intergenicRegionFeaturesFile[,1]), ncol = length(mononucleotide) + 
														length(dinucleotide) + 1);
		  
		  # --------------------------------- Features vector generate for intergenicRegion set start here ------------------------ 
		  
		  for(i in 1 : length(intergenicRegionFeaturesFile[,1]))
		  {
			tempMatrixMonoDiNucleotideComposition[i,1] = toString("0");
			
			# Use to count temp cloumn number
			tempColumnCount = 2;
			
			# -------------------------------------- mono-nucleotide assign here for intergenicRegion set -------------------------------------
			for(j in 1 : length(mononucleotide))
			{
			  columnNoMono = match(mononucleotide[j],names(intergenicRegionFeaturesFile));
			  tempMatrixMonoDiNucleotideComposition[i,tempColumnCount] = paste0(toString(tempColumnCount - 1),":", toString(intergenicRegionFeaturesFile[i,columnNoMono]));
			  tempColumnCount = tempColumnCount + 1;
			}
			
			# -------------------------------------- di-nucleotide assign here for intergenicRegion set---------------------------------------
			for(j in 1 : length(dinucleotide))
			{
			  columnNoDi = match(dinucleotide[j],names(intergenicRegionFeaturesFile));
			  tempMatrixMonoDiNucleotideComposition[i,tempColumnCount] = paste0(toString(tempColumnCount - 1),":", toString(intergenicRegionFeaturesFile[i,columnNoDi]));
			  tempColumnCount = tempColumnCount + 1;
			}
			
		  }
		  
		  # --------------------------------- Features vector generate for intergenicRegion set end here !!!------------------------ 
		  
		  tempMatrixMonoDiNucleotideComposition
		  outputFile <- paste0(toString(outputFileNameFromUser),"MonoDINucleotideCompositionSVMData",".txt");
		  write.table(as.data.frame(tempMatrixMonoDiNucleotideComposition),file = outputFile,sep=" ",quote = FALSE,col.names=FALSE, row.names = FALSE);
	  }else if(isMonoNucleotide && !isDiNucleotide && !isTriNucleotide){
		
		# This section will generate mono nucleotide composition (4) features along with indicator of class (+1 or -1)
		tempMatrixMonoNucleotideComposition = matrix(0, nrow = length(intergenicRegionFeaturesFile[,1]), ncol = length(mononucleotide) + 1);
		
		# --------------------------------- Features vector generate for intergenicRegion set start here ------------------------ 
		
		for(i in 1 : length(intergenicRegionFeaturesFile[,1]))
		{
		  tempMatrixMonoNucleotideComposition[i,1] = toString("0");
		  
		  # Use to count temp cloumn number
		  tempColumnCount = 2;
		  
		  # -------------------------------------- mono-nucleotide assign here for intergenicRegion set -------------------------------------
		  for(j in 1 : length(mononucleotide))
		  {
			columnNoMono = match(mononucleotide[j],names(intergenicRegionFeaturesFile));
			tempMatrixMonoNucleotideComposition[i,tempColumnCount] = paste0(toString(tempColumnCount - 1),":", toString(intergenicRegionFeaturesFile[i,columnNoMono]));
			tempColumnCount = tempColumnCount + 1;
		  }
		}
		
		# --------------------------------- Features vector generate for intergenicRegion set end here !!!------------------------
		
		
		tempMatrixMonoNucleotideComposition
		outputFile <- paste0(toString(outputFileNameFromUser),"MonoNucleotideCompositionSVMData", ".txt");
		write.table(as.data.frame(tempMatrixMonoNucleotideComposition),file = outputFile,sep=" ",quote = FALSE,col.names=FALSE, row.names = FALSE);
		  
	  }else if(isDiNucleotide && !isMonoNucleotide && !isTriNucleotide){
		
		# This section will generate di nucleotide composition (16) features along with indicator of class (+1 or -1)
		tempMatrixDiNucleotideComposition = matrix(0, nrow = length(intergenicRegionFeaturesFile[,1]), ncol = length(dinucleotide) + 1);
		
		# --------------------------------- Features vector generate for intergenicRegion set start here ------------------------ 
		for(i in 1 : length(intergenicRegionFeaturesFile[,1]))
		{
		  tempMatrixDiNucleotideComposition[i,1] = toString("0");
		  
		  # Use to count temp cloumn number
		  tempColumnCount = 2;
		  
		  # -------------------------------------- di-nucleotide assign here for intergenicRegion set---------------------------------------
		  for(j in 1 : length(dinucleotide))
		  {
			columnNoDi = match(dinucleotide[j],names(intergenicRegionFeaturesFile));
			tempMatrixDiNucleotideComposition[i,tempColumnCount] = paste0(toString(tempColumnCount - 1),":", toString(intergenicRegionFeaturesFile[i,columnNoDi]));
			tempColumnCount = tempColumnCount + 1;
		  }
		}
		
		# --------------------------------- Features vector generate for intergenicRegion set end here !!!------------------------ 
		
		tempMatrixDiNucleotideComposition
		outputFile <- paste0(toString(outputFileNameFromUser),"DiNucleotideCompositionSVMData", ".txt");
		write.table(as.data.frame(tempMatrixDiNucleotideComposition),file = outputFile,sep=" ",quote = FALSE,col.names=FALSE, row.names = FALSE);
		  
	  }else if(isTriNucleotide && !isMonoNucleotide && !isDiNucleotide){
		
		# This section will generate tri nucleotide composition (64) features along with indicator of class (+1 or -1)
		tempMatrixTriNucleotideComposition = matrix(0, nrow = length(intergenicRegionFeaturesFile[,1]), ncol = length(trinucleotide) + 1);
		
		
		# --------------------------------- Features vector generate for intergenicRegion set start here ------------------------ 
		for(i in 1 : length(intergenicRegionFeaturesFile[,1]))
		{
		  tempMatrixTriNucleotideComposition[i,1] = toString("0");
		  
		  # Use to count temp cloumn number
		  tempColumnCount = 2;
		  # -------------------------------------- tri-nucleotide assign here for intergenicRegion set ---------------------------------------
		  for(j in 1 : length(trinucleotide))
		  {
			columnNoTri = match(trinucleotide[j],names(intergenicRegionFeaturesFile));
			tempMatrixTriNucleotideComposition[i,tempColumnCount] = paste0(toString(tempColumnCount - 1),":", toString(intergenicRegionFeaturesFile[i,columnNoTri]));
			tempColumnCount = tempColumnCount + 1;
		  } 
		}
		# --------------------------------- Features vector generate for intergenicRegion set end here !!!------------------------ 
		tempMatrixTriNucleotideComposition
		outputFile <- paste0(toString(outputFileNameFromUser),"TriNucleotideCompositionSVMData",".txt");
		write.table(as.data.frame(tempMatrixTriNucleotideComposition),file = outputFile,sep=" ",quote = FALSE,col.names=FALSE, row.names = FALSE);
		  
	  }else{
		print("Features vectors were not generated in this case!!!!!!!!!!!!!!!!!");
	  }
  }


