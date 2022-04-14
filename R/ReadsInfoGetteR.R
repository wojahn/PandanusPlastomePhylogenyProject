ReadsInfoGetteR <- function(TN,IntermediateDirectoryPath,EstPlastLength,Trimmed)
{
  if(Trimmed == F)
  {
    OUT <- as.data.frame(matrix(nrow=length(TN),ncol=5))
    colnames(OUT) <- c("Taxon","NumberReads","InsertSize","Gbp","EstimatedCoverage")
    OUT[,1] <- TN
    ID_dirs <- list.dirs(IntermediateDirectoryPath)
    for(i in 1:length(TN))
    {
      if(length(TN) > 1)
      {
        message(sprintf("Processing sample No. %s of %s",i,length(TN)))
      }
      message("Finding files")
      this_dir <- ID_dirs[grepl(TN[i],ID_dirs)]
      this_dir_files <- list.files(this_dir)
      Unzippd <- this_dir_files[!grepl("fq.gz",this_dir_files)]
      Unzippd <- Unzippd[!grepl("Interleaved",Unzippd)]
      First_Unzippd <- Unzippd[grepl("R1",Unzippd)]
      if(length(First_Unzippd) == 0)
      {
        stop("File missing!")
      }
      message("Calculating number of reads, may take up to 1 minute!")
      numL <- as.numeric(system(sprintf("wc -l < %s/%s",this_dir,First_Unzippd), intern = T))
      message("Calculating insert length")
      headF <- system(sprintf("head -n 3 %s/%s",this_dir,First_Unzippd), intern = T)
      headF <- headF[2]
      lengthI <- nchar(headF)
      numGBP <- 2*(numL/4*lengthI)/1000000000
      numR <- numL/4*2
      coverageE <- 2*(numL/4*lengthI)/EstPlastLength
      OUT[i,2] <- numR
      OUT[i,3] <- lengthI
      OUT[i,4] <- numGBP
      OUT[i,5] <- sprintf("%s X",coverageE)
    }
  }else if(Trimmed == T){
    OUT <- as.data.frame(matrix(nrow=length(TN),ncol=5))
    colnames(OUT) <- c("Taxon","NumberReads","InsertSize","Gbp","EstimatedCoverage")
    OUT[,1] <- TN
    ID_dirs <- list.dirs(IntermediateDirectoryPath)
    for(i in 1:length(TN))
    {
      if(length(TN) > 1)
      {
        message(sprintf("Processing sample No. %s of %s",i,length(TN)))
      }
      message("Finding files")
      this_dir <- ID_dirs[grepl(TN[i],ID_dirs)]
      this_dir_files <- list.files(this_dir)
      First_Unzippd <- as.vector(this_dir_files[grepl("Interleaved",this_dir_files)])
      First_Unzippd <- First_Unzippd[grepl("fastq",First_Unzippd)]
      if(length(First_Unzippd) == 0)
      {
        stop("File missing!")
      }
      message("Calculating number of reads, may take up to 1 minute!")
      numL <- as.numeric(system(sprintf("wc -l < %s/%s",this_dir,First_Unzippd), intern = T))
      message("Calculating insert length")
      headF <- system(sprintf("head -n 3 %s/%s",this_dir,First_Unzippd), intern = T)
      headF <- headF[2]
      lengthI <- nchar(headF)
      numGBP <- (numL/4*lengthI)/1000000000
      numR <- numL/4
      coverageE <- (numL/4*lengthI)/EstPlastLength
      OUT[i,2] <- numR
      OUT[i,3] <- lengthI
      OUT[i,4] <- numGBP
      OUT[i,5] <- sprintf("%s X",coverageE)
    }
  }else{
    stop("Trimmed must be Boolean!")
  }
  return(OUT)
}
