#' Make ML/Bayesian CoPhyloPlot Trees
#' @param ML_path character: path to maximum likelihood tree file
#' @param Bayes_path character: path to bayesian inference tree file
#' @param Key_path character: path to general sample information key
#' @export
# Written by John M. A. Wojahn December 2021
# This is Free and Open-Source Software (F.O.S.S.)
# © J.M.A. Wojahn et al.
# Provided under the GNU AGPLv3 License
# Funded by EPSCoR GEM3 at Boise State University
ML_Bayes_Cophylo_MakeR <- function(ML_path,Bayes_path,Key_path)
{
  # Read in tree (need to use treeio to get metadata)
  B_Tree <- treeio::read.mrbayes(Bayes_path)
  # Convert into phylo object
  PB <- treeio::as.phylo(B_Tree)
  # Read in key
  samplesinfo <- read.csv(Key_path)
  # Read in tree
  ML_Tree <- ape::read.tree(ML_path)
  # Order nodes
  PB_Tree <- ape::ladderize(PB, right = T)
  # Order nodes
  ML_Tree <- ape::ladderize(ML_Tree, right = T)
  # Simplify labels
  PB_splitz <- strsplit(PB_Tree$tip.label,split="-")
  PBnamez <- rep(NA,length(PB_splitz))
  for(i in 1:length(PB_splitz))
  {
    if(PB_splitz[[i]][2] == "Pandanus_madagascariensis")
    {
      PBnamez[i] <- paste0(PB_splitz[[i]][2],"_",PB_splitz[[i]][1])
    }else{
      PBnamez[i] <- PB_splitz[[i]][2]
    }
  }
  PB_Tree$tip.label <- PBnamez
  Csplitz <- strsplit(ML_Tree$tip.label,split="-")
  Cnamez <- rep(NA,length(Csplitz))
  for(i in 1:length(Csplitz))
  {
    if(Csplitz[[i]][2] == "Pandanus_madagascariensis")
    {
      Cnamez[i] <- paste0(Csplitz[[i]][2],"_",Csplitz[[i]][1])
    }else{
      Cnamez[i] <- Csplitz[[i]][2]
    }
  }
  ML_Tree$tip.label <- Cnamez
  # Construct nodal association matrix
  asso <- as.matrix(data.frame(sort(ML_Tree$tip.label),sort(PB_Tree$tip.label)))
  colnames(asso) <- NULL
  Tpz <- ML_Tree$tip.label
  Tpz[Tpz == "Pandanus_sermolliana"] <- "Pandanus_sermollianus"
  ML_Tree$tip.label <- Tpz

  Tpz <- PB_Tree$tip.label
  Tpz[Tpz == "Pandanus_sermolliana"] <- "Pandanus_sermollianus"
  PB_Tree$tip.label <- Tpz
  asso[asso == "Pandanus_sermolliana"] <- "Pandanus_sermollianus"
  # Create cophyloplot
  {
  png("Bayesian_ML_Coplot.png",width = 1000, height = 1000)
  tree <- phytools::cophylo(ML_Tree,PB_Tree,assoc = asso,rotate = T)
  plot(tree)
  dev.off()
  }
}
