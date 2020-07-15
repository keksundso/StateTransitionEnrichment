#script to map cromHMM states to genes.
#todo
# biomart in function packen option einfügen das man es nicht neu laden muss

#in chromHMM they are called [CELLTYP]_[STATE-NR]_segments.bed
library(readr)
library(biomaRt)
#library(dplyr)
library(data.table)
options(scipen = 999)
#Gene Infos

#ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
#annot_df<-getBM(c("ensembl_gene_id", "external_gene_name", "chromosome_name", "start_position", "end_position","gene_biotype"), mart=ensembl)
#save(annot_df, file = "/media/ubuntuDisk/home/knalli/Dropbox/Phd/PhD-Work_Drop/h_vs_p_rnaseq_chromhmm_2020/StateTransitions/transitionEnrichment/ZS/annot_df_20200512.Rda")
load("/media/ubuntuDisk/home/knalli/Dropbox/Phd/PhD-Work_Drop/h_vs_p_rnaseq_chromhmm_2020/StateTransitions/transitionEnrichment/ZS/annot_df_20200512.Rda")

annot_df <- annot_df[annot_df$gene_biotype == "protein_coding" , c("ensembl_gene_id", "external_gene_name", "chromosome_name", "start_position", "end_position")]
annot_df$chr.annot <- paste0("chr",annot_df$chromosome_name)
annot_df$chromosome_name <- NULL
names(annot_df)[names(annot_df) == 'start_position'] <- 'start.annot'
names(annot_df)[names(annot_df) == 'end_position'] <- 'end.annot'

annot_df$ID.annot <- paste0(annot_df$chr.annot,"-",annot_df$start.annot,"-",annot_df$end.annot)

annot_bed <- annot_df[,c("chr.annot","start.annot","end.annot")]

readStateBed <-function(stateBedPath, relTime) {

  
  state_df <- read_delim(stateBedPath, 
                         "\t", escape_double = FALSE, col_names = c("chr","start","end","state"), 
                         trim_ws = TRUE)
  state_df <- state_df[(state_df$start != 0 ), ] # remove beginning of chromHMM since its urealaible
  state_df$ID <- paste0(state_df$chr,"-",state_df$start,"-",state_df$end)
  colnames(state_df) <- paste0(colnames(state_df),".",relTime)
  state_bed <- state_df[,1:3]
  
  return(list(state_df,state_bed))
}

mergeStateFiles <- function(stateListA, stateListB){
  df_A <- as.data.table(stateListA[[2]])
  df_B <- as.data.table(stateListB[[2]])
  
  setkeyv(df_A, names(df_A))
  setkeyv(df_B, names(df_B))
  ans <- foverlaps(df_B, df_A, nomatch=0L)#,minoverlap=10) #is not implimented yet, will do it by below

  df <- as.data.frame((ans))
  
  df[c('start.overlap', 'end.overlap')] <- t(apply(df, 1, function(x) 
    sort(x)[c(length(x)-3, length(x) - 2)]))
  
  # minoverlap by hand
  df$start.overlap <- as.integer(as.character(df$start.overlap))
  df$end.overlap <- as.integer(as.character(df$end.overlap))
  df <- df[abs(df$start.overlap - df$end.overlap)> 10 , ]
  
  df$ID.A <- paste0(df$chr.B,"-",df$start.A,"-",df$end.A)
  df$ID.B <- paste0(df$chr.B,"-",df$start.B,"-",df$end.B)
  
  df <- merge(df, stateListA[[1]][,c("state.A","ID.A")], by ="ID.A")
  df <- merge(df, stateListB[[1]][,c("state.B","ID.B")], by ="ID.B")
  df$transiton <- paste(substring(df$state.A,2),substring(df$state.B,2) ,sep = ".")

  transiton_df <-  df[,c("chr.B","start.overlap","end.overlap","transiton")]
  names(transiton_df)[names(transiton_df) == 'chr.B'] <- 'chr.overlap'  
  transiton_df$ID.Overlap <- paste0(transiton_df$chr.overlap,"-",transiton_df$start.overlap,"-",transiton_df$end.overlap)
  transiton_bed <- transiton_df[, c("chr.overlap","start.overlap","end.overlap")]
  
  return(list(transiton_df, transiton_bed))
}  

creatVisTransbed <- function(transiton_df, outPutPath = getwd()){

  if("external_gene_name"%in% colnames(transiton_df)){
         tmp <- transiton_df[,c("chr.overlap","start.overlap","end.overlap","transiton","external_gene_name")]
         outPutPath <- file.path(outPutPath,"transitionGeneVis.bed")
         }
  if(!("external_gene_name"%in% colnames(transiton_df))){
       tmp <- transiton_df[,c("chr.overlap","start.overlap","end.overlap","transiton")]
       outPutPath <- file.path(outPutPath,"transitonVis.bed")}

  rgbValue <- c()
  rbgTable <- as.data.frame(col2rgb(sample(grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)], length(unique(tmp$transiton ) ))))
  for (variable in colnames(rbgTable)) {
    rgbValue <-c(rgbValue, paste(rbgTable[,variable][1],rbgTable[,variable][2],rbgTable[,variable][3], sep = ","))
  }
  tmpA <- data.frame(transiton = unique(tmp$transiton), color = rgbValue)
  

  
  tmp$ka <- 1000
  tmp$strand <- "."
  tmp$startAgain <- tmp$start.overlap
  tmp$endAgian <- tmp$end.overlap
  tmpB <- colnames(tmp)
  tmp <- merge(tmp, tmpA, by="transiton")
  tmp <- tmp[,c(tmpB, "color")]
  if("external_gene_name"%in% colnames(transiton_df)){
    tmp$transiton <- paste0(tmp$transiton,"_",tmp$external_gene_name)
    tmp$external_gene_name <- NULL
  }
  
  write.table(tmp,file =outPutPath, sep = "\t",row.names = FALSE, col.names = FALSE, quote = FALSE)
 
}

MergeTransitionWithGeneAnno <- function(annot_bed,mergedStateList){
  df_A <- as.data.table(annot_bed)
  df_B <- as.data.table(mergedStateList[[2]])
  
  setkeyv(df_A, names(df_A))
  setkeyv(df_B, names(df_B))
  ans <- foverlaps(df_B, df_A, nomatch=0L)
  
  ans$ID.annot <- paste0(ans$chr.overlap,"-",ans$start.annot,"-",ans$end.annot)
  ans$ID.Overlap <- paste0(ans$chr.overlap,"-",ans$start.overlap,"-",ans$end.overlap)
  ans <- ans[,c("ID.annot","ID.Overlap")]
  ans <- merge(annot_df, ans, by="ID.annot")
  ans <- merge(ans, mergedStateList[[1]], by="ID.Overlap")
  creatVisTransbed(ans)
 

  return(ans)
}


stateListA <- readStateBed("StateTransitions/chromatinStates/chromHMM/outputData/myLernOutPutBig/PC_8_segments.bed","A")
stateListB <- readStateBed("StateTransitions/chromatinStates/chromHMM/outputData/myLernOutPutBig/HC_8_segments.bed","B")
mergedStateList <- mergeStateFiles(stateListA,stateListB)
creatVisTransbed(mergedStateList[[1]])
transitionGeneDF <- MergeTransitionWithGeneAnno(annot_bed,mergedStateList)
