DeIUPACifieR <- function(file)
{
  badz <- c(letters,LETTERS,"?")
  goodz <- c("A","T","G","C","a","t","g","c")
  badz <- badz[!badz %in% goodz]
  linez <- readLines(file)
  linez_fixed <- linez
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                     max = length(linez), # Maximum value of the progress bar
                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
                     width = 50,   # Progress bar width. Defaults to getOption("width")
                     char = "=")   # Character used to create the bar
  for(i in 1:length(linez))
  {
    setTxtProgressBar(pb, i)
    if(!grepl(">",linez[i]))
    {
      splitz <- unlist(strsplit(linez[i], split = ""))
      if(T %in% (badz %in% splitz))
      {
        splitz[splitz %in% badz] <- "N"
        gluedz <- paste(splitz, collapse="")
        linez_fixed[i] <- gluedz
      }
    }
  }
  close(pb)
  filenude <- gsub(".fasta","",file)
  write.table(linez_fixed,sprintf("%s_No_IUPAC.fasta",filenude),row.names = F, col.names = F, quote = F)
}
