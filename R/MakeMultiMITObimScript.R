MakeMultiMITObimScript <- function(FullMITObimPath, FullMIRAPath, FullSeedPath, TN, ScaffoldsPath,removeits,parallel)
{
  filez <- list.files()
  if("MasterMITObimScript.sh" %in% filez)
  {
    message("Deleting preëxisting MasterMITObimScript")
    file.remove("MasterMITObimScript.sh")
  }
  for(i in 1:length(TN))
  {
    message(sprintf("Processing No. %s of %s",i,length(TN)))
    message("Making MITObim output directory (MITObim_Genomes)")
    syscom <- sprintf("mkdir -p MITObim_Genomes/%s", TN[i])
    system(syscom)
    message("Listing files")
    filez <- list.files(sprintf("Intermediate_Data/%s", TN[i]))
    message("Identifying interleaved fasta")
    combo_filez <- filez[grep("Interleaved",filez)]
    combo_filez <- combo_filez[grep(".fasta",combo_filez)]
    MBfilez <- list.files(sprintf("MITObim_Genomes/%s", TN[i]))
    if(T %in% grepl("Interleaved",MBfilez))
    {
      message("Interleaved file already there...moving on!")
    }else{
      message("Moving interleaved fastq to MITObim_Genomes")
      filez <- list.files(sprintf("Intermediate_Data/%s", TN[i]))
      combo_filez <- filez[grep("Interleaved",filez)]
      combo_filez <- combo_filez[grep(".fastq",combo_filez)]
      syscom <- sprintf("cp Intermediate_Data/%s/%s MITObim_Genomes/%s/%s",TN[i],combo_filez,TN[i],combo_filez)
      system(syscom)
    }
    if(T %in% grepl("ThisSeed.fasta",MBfilez))
    {
      message("Reference file already there...moving on!")
    }else{
      message("Copying reference to MITObim_Genomes")
      syscom <- sprintf("cp %s MITObim_Genomes/%s/ThisSeed.fasta",FullSeedPath[i],TN[i])
      system(syscom)
    }
    message("Compiling and writing MITObim script")
    message("Please ignore any warnings :)")
    system(sprintf("rm MITObim_Genomes/%s/MITObimScript.sh",TN[i]))
    commande <- sprintf("cd MITObim_Genomes/%s",TN[i])
    syscom <- sprintf("echo '%s' >> MITObim_Genomes/%s/MITObimScript.sh",commande,TN[i])
    system(syscom)
    commande <- sprintf("%s -start 0 -end 30 -sample %s -ref genome -readpool %s_Interleaved.fastq --quick ThisSeed.fasta -mirapath %s --clean --verbose --paired",FullMITObimPath,TN[i],TN[i],FullMIRAPath)
    syscom <- sprintf("echo '%s' >> MITObim_Genomes/%s/MITObimScript.sh",commande,TN[i])
    system(syscom)
    message("Setting permissions")
    syscom <- sprintf("chmod u+x MITObim_Genomes/%s/MITObimScript.sh",TN[i])
    system(syscom)
    message("Finishing script entry")
    syscom <- sprintf("./MITObim_Genomes/%s/MITObimScript.sh",TN[i])
    commande <- sprintf("echo '%s' >> MasterMITObimScript.sh",syscom)
    system(commande)
    if(removeits == T)
    {
      message("Cleaning folder of old iterations")
      if(T %in% grepl("iteration",MBfilez))
      {
        toremove <- MBfilez[grepl("iteration",MBfilez)]
        for(j in 1:length(toremove))
        {
          commande <- sprintf("rm -rf MITObim_Genomes/%s/%s",TN[i],toremove[j])
          system(commande)
        }
      }
    }else{
      message("Keeping old iterations, as requested")
    }
  }
  syscom <- sprintf("chmod u+x MasterMITObimScript.sh")
  system(syscom)
  if(!is.null(parallel))
  {
    MMbS <- readLines("MasterMITObimScript.sh")
    num <- ceiling(length(MMbS)/parallel)
    nomenpar <- c()
    countR <- 0
    for(i in 1:parallel)
    {
      if(i == 1)
      {
        them <- MMbS[seq(from = countR+1, to = countR+num, by = 1)]
        countR <- num
        write.table(them, file = sprintf("MasterMITObimScriptPar_%s.sh",i),row.names = F,sep="\t", quote = FALSE,col.names = F)
        nomenpar <- c(nomenpar,sprintf("MasterMITObimScriptPar_%s.sh &",i))
      }else if(i == parallel){
        them <- MMbS[seq(from = countR+1, to = length(MMbS), by = 1)]
        write.table(them, file = sprintf("MasterMITObimScriptPar_%s.sh",i),row.names = F,sep="\t", quote = FALSE,col.names = F)
        nomenpar <- c(nomenpar,sprintf("MasterMITObimScriptPar_%s.sh",i))
      }else{
        them <- MMbS[seq(from = countR+1, to = countR+num, by = 1)]
        countR <- countR+num
        write.table(them, file = sprintf("MasterMITObimScriptPar_%s.sh",i),row.names = F,sep="\t", quote = FALSE,col.names = F)
        nomenpar <- c(nomenpar,sprintf("MasterMITObimScriptPar_%s.sh &",i))
        syscom <- sprintf("chmod u+x MasterMITObimScriptPar_%s.sh",i)
        system(syscom)
      }
    }
    write.table(as.character(sprintf("./%s> log%s %s",nomenpar,seq(from = 1, to = parallel),c(rep("&",parallel-1),""))), file = sprintf("MasterParallelMITObimScript.sh"),row.names = F,sep="\t", quote = FALSE,col.names = F)
    syscom <- sprintf("chmod u+x MasterParallelMITObimScript.sh")
    system(syscom)
    message("All done")
    message("Please run MasterParallelMITObimScript.sh NOT the MITObimScript.sh in the subdirectory(ies)!")
  }else{
    message("All done")
    message("Please run MasterMITObimScript.sh NOT the MITObimScript.sh in the subdirectory(ies)!")
  }
}
