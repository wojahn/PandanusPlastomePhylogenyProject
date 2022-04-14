MITObimStatsGetteR <- function(tolook)
{
  message("Listing directories")
  taxonz <- list.files(tolook)
  exceptionz <-  c("Pandan0016","Pandan0018","Pandan0022",
                   "Pandan0030","Pandan0042","Pandan0049",
                   "Pandan0056")
  for(i in 1:length(exceptionz))
  {
    taxonz <- taxonz[!grepl(exceptionz[i],taxonz)]
  }
  OUT <- as.data.frame(matrix(nrow=length(taxonz), ncol=6))
  OUT[,1] <- taxonz
  colnames(OUT) <- c("Taxon","Length","AvgCoverage","NumReads","NumberNs","GC")
  for(i in 1:length(taxonz))
  {
    message(sprintf("Processing directory no. %s of %s",i,length(taxonz)))
    message("Listing files")
    filez <- list.files(sprintf("%s/%s",tolook,taxonz[i]))
    itz <- filez[grepl("iteration",filez)]
    message("Entering most recent iteration")
    Numeroitz <- as.numeric(gsub("iteration","",itz))
    sortedz <- itz[which(Numeroitz == max(Numeroitz))]
    message("Checking if logfile exists")
    if(file.exists(sprintf("%s/%s/%s/%s-genome_assembly/%s-genome_d_info/%s-genome_info_contigstats.txt",tolook,taxonz[i],sortedz,taxonz[i],taxonz[i],taxonz[i])))
    {
      message("Reading logfile")
      statz <- read.delim(sprintf("%s/%s/%s/%s-genome_assembly/%s-genome_d_info/%s-genome_info_contigstats.txt",tolook,taxonz[i],sortedz,taxonz[i],taxonz[i],taxonz[i]), sep = "\t")
      message("Processing logfile")
      message("Writing info")
      OUT[i,2] <- statz[1,colnames(statz) == "length"]
      OUT[i,3] <- statz[1,colnames(statz) == "av.cov"]
      OUT[i,4] <- statz[1,colnames(statz) == "X..reads"]
      OUT[i,5] <- statz[1,colnames(statz) == "CnN"]
      OUT[i,6] <- statz[1,colnames(statz) == "GC."]
    }else{
      message("No logfile, MITObim failed, noting in OUT")
      OUT[i,2] <- "None"
      OUT[i,3] <- "None"
      OUT[i,4] <- "None"
      OUT[i,5] <- "None"
      OUT[i,6] <- "None"
    }
  }
  message("Done!")
  return(OUT)
}
