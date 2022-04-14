#' Make ML/Bayesian Combination Tree With Labels
#' @param ML_path character: path to maximum likelihood tree file
#' @param Bayes_path character: path to bayesian inference tree file
#' @param ISO_path character: path to ISO country codes key
#' @param Key_path character: path to general sample information key
#' @param Subgen_path character: path to subgeneric taxonomy key
#' @export
# Written by John M. A. Wojahn December 2021
# This is Free and Open-Source Software (F.O.S.S.)
# © J.M.A. Wojahn et al.
# Provided under the GNU AGPLv3 License
# Funded by EPSCoR GEM3 at Boise State University
ML_Bayes_Combo_Tree_MakeR <- function(ML_path,Bayes_path,ISO_path,Key_path,Subgen_path)
{
  # Read in tree
  ML_Tree <- ape::read.tree(ML_path)
  # Order nodes
  ML_Tree <- ape::ladderize(ML_Tree, right = T)
  # Simplify labels
  laboz <- ML_Tree$tip.label
  labozsplitz <- strsplit(laboz, split = "-")
  labozgoodz <- rep(NA,length(labozsplitz))
  for(i in 1:length(labozsplitz))
  {
   labozgoodz[i] <- labozsplitz[[i]][2]
  }
  labozgoodz <- gsub("_"," ",labozgoodz)
  ML_Tree$tip.label <- labozgoodz
 # ML_Tree$node.label
 # ML_Tree$edge
  # Read in tree (need to use treeio to get metadata)
  require(treeio)
  B_Tree <- treeio::read.mrbayes(Bayes_path)
  # Convert into phylo object
  PB <- treeio::as.phylo(B_Tree)
  # Order nodes
  PB <- ape::ladderize(PB, right = T)
  # Simplify labels
  laboz <- PB$tip.label
  labozsplitz <- strsplit(laboz, split = "-")
  labozgoodz <- rep(NA,length(labozsplitz))
  for(i in 1:length(labozsplitz))
  {
   labozgoodz[i] <- labozsplitz[[i]][2]
  }
  labozgoodz <- gsub("_"," ",labozgoodz)
  PB$tip.label <- labozgoodz
  # Extract posterior probabilities
  PostProbz <- B_Tree@data[["prob"]]
  message(sprintf("Mean BPP: %s",mean(as.numeric(PostProbz))))
  PostProbzNodez <- B_Tree@data[["node"]]
  tipz <- ML_Tree$tip.label
  nodelabz <- ML_Tree$node.label
  message(sprintf("Mean BS: %s",mean(as.numeric(nodelabz[nodelabz != ""]))))
  combonodez <- rep(NA,ML_Tree$Nnode)
  # Combine posterior probabilities with bootstrap percentages using MRCAs of nodes
  for(i in 1:(length(tipz)-1))
  {
    MLnode <- phytools::findMRCA(ML_Tree,c(tipz[i],tipz[i+1]),type="node")
    ML_boot <- ML_Tree$node.label[MLnode - Ntip(ML_Tree)]
    BayesNode <- phytools::findMRCA(PB,c(tipz[i],tipz[i+1]),type="node")
    Bayes_PP <- round(as.numeric(PostProbz[PostProbzNodez == BayesNode]),2)
    if(length(Bayes_PP) == 0)
    {
      Bayes_PP <- "NA"
    }
    if(ML_boot == "")
    {
      ML_boot <- "NA"
    }
    combonodez[MLnode - Ntip(ML_Tree)] <- sprintf("%s/%s", ML_boot, Bayes_PP)
  }
  ComboTree <- ML_Tree
  ComboTree$node.label <- combonodez

  # Make labels
  iso_countries <- read.csv(ISO_path)
  samplesinfo <- read.csv(Key_path)
  subgencodes <- read.csv(Subgen_path)
  rawlabz <- ComboTree$tip.label
  rawlabz_big <- rep(NA,length(rawlabz))
  for(i in 1:length(rawlabz))
  {
    littleinfo <- samplesinfo[samplesinfo$Taxon ==
                                gsub(" ","_",rawlabz[i]),]
    littleinfosmol <- unique(littleinfo[,c(15,16,20)])
    subgenabbr <- subgencodes[subgencodes$Subgenus ==
                                littleinfosmol[1,1],2]
    if(length(subgenabbr)==0)
    {
      subgenabbr <- "Unknown"
    }
    sectionfull <- littleinfosmol[,2]
    country <- littleinfosmol[,3]
    isocountry <- iso_countries[iso_countries$Country == country,3]
    if(length(isocountry) == 0)
    {
      isocountry <- iso_countries[gsub(" .*","",iso_countries$Country) == country,3]
      if(length(isocountry) == 0)
      {
        if(country == "Reunion" || country == "Réunion")
        {
          isocountry <- "REU"
        }else if(country == "Borneo"){
          isocountry <- "MYS"
        }else if(country == "New Guinea"){
          isocountry <- "PNG"
        }else{
          isocountry <- country
        }
      }
    }
    rawlabz_big[i] <- sprintf("%s_%s_%s_%s",gsub(" ","_",rawlabz[i]),subgenabbr,sectionfull,isocountry)
  }
  ComboTree$tip.label <- rawlabz_big

  ComboTreeFull <- ComboTree
  tipz <- ComboTree$tip.label
  splitz <- strsplit(tipz,split="_")
  smalltipz <- tipz[sapply(splitz, "[[", 1) %in% c("Sararanga","Freycinetia","Benstonea")]

  ComboTree <- ape::drop.tip(ComboTree,smalltipz)
  # Plot combination tree


  SCII_node <- castor::get_pairwise_mrcas(ComboTree, "Pandanus_multispicatus_VIN_Ignotus_SYC",
                             "Pandanus_aquaticus_PAN_Semikeura_AUS" ,
                             check_input=TRUE)

  SCI_node <- castor::get_pairwise_mrcas(ComboTree,
                                "Pandanus_irregularis_RYK_Mydiophylla_THA",
                             "Pandanus_zea_LOP_Maysops_AUS",
                             check_input=TRUE)
  Tpz <- ComboTree$tip.label
  Tpz[Tpz == "Pandanus_sermolliana_PAN_Pandanus_MDG"] <- "Pandanus_sermollianus_PAN_Pandanus_MDG"
  ComboTree$tip.label <- Tpz
  ComboTree <- ape::ladderize(ComboTree)
  {
    pdf("ComboTreeUpdated.pdf",width=20,height=15)
    ape::plot.phylo(ComboTree, cex = 0.7,label.offset=0.00009,x.lim	=0.015)
    ape::nodelabels(as.character(ComboTree[["node.label"]]), cex = 0.55,frame="none",bg="white",adj=c(0,0),col="blue")
    #phytools::cladelabels(ComboTree,"Clade II",SCII_node,cex=1,offset=-4,wing.length=0,orientation="vertical")
    #phytools::cladelabels(ComboTree,"Clade I",SCI_node,cex=1,offset=-0.60,wing.length=0,orientation="vertical")
    #ape::add.scale.bar()
    dev.off()
  }
}
