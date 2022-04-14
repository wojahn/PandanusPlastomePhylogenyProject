Do_Multi_Trimmomatic_Command <- function(TN)
{
  for(i in 1:length(TN))
  {
    print(sprintf("Processing No. %s of %s",i,length(TN)))
    # Make intermediate output folder
    syscom <- sprintf("mkdir Intermediate_Data/%s",TN[i])
    system(syscom)
    # Jar file path
    jarfile <- "/usr/local/Cellar/trimmomatic/0.39/libexec/trimmomatic-0.39.jar"
    # Input forward file name
    filez <- list.files(sprintf("Raw_Data/%s", TN[i]))
    f_filez <- filez[grep("R1",filez)]
    r_filez <- filez[grep("R2",filez)]
    input_f <- sprintf("Raw_Data/%s/%s", TN[i], f_filez)
    # Input reverse file name
    input_r <- sprintf("Raw_Data/%s/%s", TN[i], r_filez)
    # Output foward unpaired file name
    FFN_base <- gsub(".fastq.gz","",input_f)
    FFN_base <- gsub(sprintf("Raw_Data/%s/", TN[i]),"",FFN_base)
    RFN_base <- gsub(".fastq.gz","",input_r)
    RFN_base <- gsub(sprintf("Raw_Data/%s/", TN[i]),"",RFN_base)
    output_f_p <- sprintf("Intermediate_Data/%s/%s_f_p.fq.gz", TN[i], FFN_base)
    # Output foward paired file name
    output_f_u <- sprintf("Intermediate_Data/%s/%s_f_u.fq.gz", TN[i], FFN_base)
    # Output reverse paired file name
    output_r_p <- sprintf("Intermediate_Data/%s/%s_f_p.fq.gz", TN[i], RFN_base)
    # Output reverse unpaired file name
    output_r_u <- sprintf("Intermediate_Data/%s/%s_f_u.fq.gz", TN[i], RFN_base)
    syscom <- sprintf("java -jar %s PE -threads 14 -phred33 %s %s %s %s %s %s HEADCROP:20 SLIDINGWINDOW:4:30 MINLEN:30", jarfile, input_f, input_r, output_f_p, output_f_u, output_r_p, output_r_u)
    system(syscom)
  }
  return(0)
}
