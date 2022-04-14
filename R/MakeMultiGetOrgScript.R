MakeMultiGetOrgScript <- function(numcores,TN,NumRo,kmersizes,Template,ws,reftype)
{
  # Number of cores to utilize
  print("Setting cores")
  numcores <- numcores
  for(i in 1:length(TN))
  {
    print(sprintf("Processing No. %s of %s",i,length(TN)))
    print("Listing files")
    syscom <- sprintf("mkdir -p Genomes/%s", TN[i])
    system(syscom)
    # Input is trimmed compressed phred+33 (.fq.gz) file name
    filez <- list.files(sprintf("Intermediate_Data/%s", TN[i]),recursive=T)
    f_filez <- filez[grep("R1",filez)]
    f_filez <- f_filez[grep("p",f_filez)]
    f_filez <- f_filez[grep("gz",f_filez)]
    r_filez <- filez[grep("R2",filez)]
    r_filez <- r_filez[grep("p",r_filez)]
    r_filez <- r_filez[grep("gz",r_filez)]
    print("Getting file names")
    input_f <- sprintf("Intermediate_Data/%s/%s", TN[i], f_filez)
    input_r <- sprintf("Intermediate_Data/%s/%s", TN[i], r_filez)
    # Output is asm file name
    print("Creating output directory")
    output_name <- sprintf("Genomes/%s/%s_Assembly", TN[i], TN[i])
    if(i == 1)
    {
      system("echo '#!/bin/bash' - > GetOrganelleScript.sh") #shibboleth
    }
    print("Creating command")
    if(ws == "auto")
    {
      if(Template == "auto")
      {
        syscom <- sprintf("get_organelle_from_reads.py -1 %s -2 %s -o %s -t %s -R %s -k %s -F %s --overwrite", input_f, input_r, output_name, numcores, NumRo, kmersizes, reftype)
      }else{
         syscom <- sprintf("get_organelle_from_reads.py -1 %s -2 %s -o %s -t %s -R %s -k %s -F %s -s %s --overwrite", input_f, input_r, output_name, numcores, NumRo, kmersizes,reftype, Template)
      }
    }else{
      if(Template == "auto")
      {
        syscom <- sprintf("get_organelle_from_reads.py -1 %s -2 %s -o %s -t %s -R %s -k %s -F %s --overwrite -w %s", input_f, input_r, output_name, numcores, NumRo, kmersizes, reftype, ws)
      }else{
         syscom <- sprintf("get_organelle_from_reads.py -1 %s -2 %s -o %s -t %s -R %s -k %s -F %s -s %s --overwrite -w %s", input_f, input_r, output_name, numcores, NumRo, kmersizes, reftype, Template, ws)
      }
    }
    print("Writing command")
    syscom <- sprintf("echo '%s' >> Combo_GetOrganelleScript.sh",syscom)
    system(syscom)
  }
  print("Setting permissions")
  system("chmod u+x Combo_GetOrganelleScript.sh")
  system(syscom)
  return(0)
}
