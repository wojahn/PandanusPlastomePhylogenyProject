ResultFromLog <- function(TN,subdirect,dub)
{
  numscafz <- rep(NA,length(TN)) # make blank
  for(i in 1:length(TN)) # extract info
  {
    loc <- sprintf("%s/%s/%s_Assembly/get_org.log.txt",subdirect,TN[i],TN[i])
    logfiletext <- read.delim(loc, sep = "\t")
    justtxt <- as.character(logfiletext[grep("Result status of embplant_pt", logfiletext[,1]),])
    numscaf <- gsub(".*Result status of embplant_pt: ","",justtxt)
    numscafz[i] <- numscaf
  }
  if(dub == T) # if wish to double
  {
    j <- 1
    i <- 1
    out <- rep(NA,2*length(TN))
    while(j <= 2*length(TN))
    {
      out[j] <- numscafz[i]
      if(j %% 2 == 0) #if j is even
      {
        i <- i + 1
      }
      j <- j + 1
    }
  }else{
    out <- numscafz
  }
  return(out)
}
