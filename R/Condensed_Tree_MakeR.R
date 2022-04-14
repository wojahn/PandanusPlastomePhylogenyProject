#' Makes Tree With Simplified Condensed Pandanus Branch
#' @param ML_path character: path to maximum likelihood tree
#' @param Bayes_path character: path to Bayesian inference tree
#' @param ISO_path character: path to ISO country codes key
#' @param Key_path character: path to general samples key
#' @param Subgen_path character: path to subgenera taxonomical key
#' @export
# Written by John M. A. Wojahn December 2021
# This is Free and Open-Source Software (F.O.S.S.)
# © J.M.A. Wojahn et al.
# Provided under the GNU AGPLv3 License
# Funded by EPSCoR GEM3 at Boise State University
Condensed_Tree_MakeR <- function(ML_path, Bayes_path, ISO_path, Key_path, Subgen_path)
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
  #ML_Tree$node.label
  #ML_Tree$edge
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
  PostProbzNodez <- B_Tree@data[["node"]]
  tipz <- ML_Tree$tip.label
  nodelabz <- ML_Tree$node.label
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
  namez <- ComboTreeFull$tip.label
  namezsmall <- namez[!gsub("_.*","",namez) %in% c("Sararanga","Freycinetia","Benstonea")]
  namezsmall <- namezsmall[!namezsmall == "Pandanus_eydouxia_EYD_Eydouxia_MUS"]
  generatree <- ape::drop.tip(ComboTreeFull,namezsmall)
  #generatree$tip.label <- sapply(strsplit(generatree$tip.label,split="_"), "[[", 1)
  gttl <- generatree$tip.label
  gttl[gttl == "Pandanus_eydouxia_EYD_Eydouxia_MUS"] <- "Pandanus"
  generatree$tip.label <- gttl
  C2_node <- castor::get_pairwise_mrcas(generatree, "Benstonea_pachyphylla_IGNO_Ignotus_MYS",
                             "Benstonea_epiphytica_IGNO_Ignotus_MYS" ,
                             check_input=TRUE)

  C3c_node <- castor::get_pairwise_mrcas(generatree, "Benstonea_inquilina_IGNO_Ignotus_MYS",
                             "Benstonea_brevistyla_IGNO_Ignotus_MYS" ,
                             check_input=TRUE)
  tipz <- generatree$tip.label
  newtipz <- tipz
  for(i in 1:length(tipz))
  {
    if(tipz[i] != "Pandanus")
    {
      splitz <- unlist(strsplit(tipz[i],split="_"))
      newtipz[i] <- paste(c(splitz[1:2],splitz[length(splitz)]),collapse="_")
    }
  }
  generatree$tip.label <- newtipz
  # PLot genera tree
  {
    pdf("GeneraTree.pdf",width=12,height=5)
      ape::plot.phylo(generatree, cex = 1.5,label.offset=0.00009)
       #phytools::cladelabels(generatree,"Clade 2",C2_node,cex=1,offset=-1.3,wing.length=0,orientation="vertical")
            #phytools::cladelabels(generatree,"Clade 3a",C3c_node,cex=1,offset=-3.75,wing.length=0,orientation="vertical")
            ape::add.scale.bar()
      dev.off()
  }
}
