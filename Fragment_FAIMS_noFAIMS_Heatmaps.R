BiocManager::install("EBImage")

library(pcaMethods)
library(ggplot2)
library(pheatmap)
library(plotly)
library(reshape2)
library(RColorBrewer)
library(plyr)
library("EBImage")

#### Upload Data ####

# data.frame of Eclipse 45 min runs using HP column. Searched with new aliqouts to confirm intensity trends 
CompiledFragments_0_4_noFAIMS <- read.csv("H:/Projects/Proteomics/Zymomona/FAIMS/DataAnalysis/SummedIntensity_Infusions/AIEIVDQALDR/CompiledFragmentResults_IW_0_4_Zymo_noFAIMS.csv", 
                                            header = TRUE, sep = ",", stringsAsFactors = FALSE)


CompiledFragments_0_4_FAIMS <- read.csv("H:/Projects/Proteomics/Zymomona/FAIMS/DataAnalysis/SummedIntensity_Infusions/AIEIVDQALDR/CompiledFragmentResults_IW_0_4_Zymo_FAIMS.csv", 
                                          header = TRUE, sep = ",", stringsAsFactors = FALSE)


CompiledFragments_0_6_noFAIMS <- read.csv("H:/Projects/Proteomics/Zymomona/FAIMS/DataAnalysis/SummedIntensity_Infusions/AIEIVDQALDR/CompiledFragmentResults_IW_0_6_Zymo_noFAIMS.csv", 
                                          header = TRUE, sep = ",", stringsAsFactors = FALSE)


CompiledFragments_0_6_FAIMS <- read.csv("H:/Projects/Proteomics/Zymomona/FAIMS/DataAnalysis/SummedIntensity_Infusions/AIEIVDQALDR/CompiledFragmentResults_IW_0_6_Zymo_FAIMS.csv", 
                                        header = TRUE, sep = ",", stringsAsFactors = FALSE)

CompiledFragments_0_8_noFAIMS <- read.csv("H:/Projects/Proteomics/Zymomona/FAIMS/DataAnalysis/SummedIntensity_Infusions/AIEIVDQALDR/CompiledFragmentResults_IW_0_8_Zymo_noFAIMS.csv", 
                                          header = TRUE, sep = ",", stringsAsFactors = FALSE)


CompiledFragments_0_8_FAIMS <- read.csv("H:/Projects/Proteomics/Zymomona/FAIMS/DataAnalysis/SummedIntensity_Infusions/AIEIVDQALDR/CompiledFragmentResults_IW_0_8_Zymo_FAIMS.csv", 
                                        header = TRUE, sep = ",", stringsAsFactors = FALSE)

CompiledFragments_1_noFAIMS <- read.csv("H:/Projects/Proteomics/Zymomona/FAIMS/DataAnalysis/SummedIntensity_Infusions/AIEIVDQALDR/CompiledFragmentResults_IW_1_Zymo_noFAIMS.csv", 
                                          header = TRUE, sep = ",", stringsAsFactors = FALSE)


CompiledFragments_1_FAIMS <- read.csv("H:/Projects/Proteomics/Zymomona/FAIMS/DataAnalysis/SummedIntensity_Infusions/AIEIVDQALDR/CompiledFragmentResults_IW_1_Zymo_FAIMS.csv", 
                                        header = TRUE, sep = ",", stringsAsFactors = FALSE)


CompiledFragments_1_3_noFAIMS <- read.csv("H:/Projects/Proteomics/Zymomona/FAIMS/DataAnalysis/SummedIntensity_Infusions/AIEIVDQALDR/CompiledFragmentResults_IW_1_3_Zymo_noFAIMS.csv", 
                                        header = TRUE, sep = ",", stringsAsFactors = FALSE)


CompiledFragments_1_3_FAIMS <- read.csv("H:/Projects/Proteomics/Zymomona/FAIMS/DataAnalysis/SummedIntensity_Infusions/AIEIVDQALDR/CompiledFragmentResults_IW_1_3_Zymo_FAIMS.csv", 
                                      header = TRUE, sep = ",", stringsAsFactors = FALSE)

CompiledFragments_2_noFAIMS <- read.csv("H:/Projects/Proteomics/Zymomona/FAIMS/DataAnalysis/SummedIntensity_Infusions/AIEIVDQALDR/CompiledFragmentResults_IW_2_Zymo_noFAIMS.csv", 
                                          header = TRUE, sep = ",", stringsAsFactors = FALSE)


CompiledFragments_2_FAIMS <- read.csv("H:/Projects/Proteomics/Zymomona/FAIMS/DataAnalysis/SummedIntensity_Infusions/AIEIVDQALDR/CompiledFragmentResults_IW_2_Zymo_FAIMS.csv", 
                                        header = TRUE, sep = ",", stringsAsFactors = FALSE)

#### Format Data ####

# x = CompiledFragments_0_4_FAIMS
meta_replicates <- function(x){
  samples <- colnames(x[,-c(1:2)])
  samples <- gsub("^.*?_R","R",samples)
  
  df <- x[,-c(1:2)]
  rownames(df) <- x$Mass.Feature
  colnames(df) <- samples
  df[df == 0] <- NA
  
  df <- df[-which(rownames(df)=="b1"),]
  df <- df[-which(rownames(df)=="y1"),]
  
  df_t <- t(df)
  
  meta <- as.data.frame(stringr::str_split_fixed(rownames(df_t),"_",4))
  colnames(meta) <- c("Resolution", "MaxInjection", "AGC", "Replicate")
  
  meta[,1] <- sub('R','', meta[,1])
  
  meta[,2] <- sub('MI', '', meta[,2])
  meta[,2] <- sub(502, 0.502, meta[,2])
  meta[,2] <- as.factor(as.numeric(meta[,2]))
  
  meta[,4] <- as.numeric(sub("*.raw",'', meta[,4]))
  
  rownames(meta) <- rownames(df_t)
  meta <- meta[order(meta[,1],meta[,2]),]
  
  return(meta)
  
}


# x = CompiledFragments_0_4_FAIMS
# m = meta_0_4_noFAIMS
df_fragments_log2 <- function(x,m){
  samples <- colnames(x[,-c(1:2)])
  samples <- gsub("^.*?_R","R",samples)
  
  df <- x[,-c(1:2)]
  rownames(df) <- x$Mass.Feature
  colnames(df) <- samples
  df[df == 0] <- NA
  
  df_log2 <- log2(df)
  df_log2[is.na(df_log2)] <- 0
  
  df_log2 <- df_log2[-which(rownames(df_log2)=="b1"),]
  df_log2 <- df_log2[-which(rownames(df_log2)=="y1"),]
  
  for(i in 1:ncol(df_log2)){
    mostIntense = max(df_log2[,i])
    
    for(j in 1:nrow(df_log2)){
      df_log2[j,i] <- df_log2[j,i]/mostIntense
    }
  }

  df_log2 <- df_log2[rownames(m)]
  
  return(df_log2)
  
}

#### Meta and Log2 dataframes ####

meta_0_4_noFAIMS <- meta_replicates(CompiledFragments_0_4_noFAIMS)
df_0_4_noFAIMS_log2 <- df_fragments_log2(CompiledFragments_0_4_noFAIMS,meta_0_4_noFAIMS)

meta_0_4_FAIMS <- meta_replicates(CompiledFragments_0_4_FAIMS)
df_0_4_FAIMS_log2 <- df_fragments_log2(CompiledFragments_0_4_FAIMS, meta_0_4_FAIMS)


meta_0_6_noFAIMS <- meta_replicates(CompiledFragments_0_6_noFAIMS)
df_0_6_noFAIMS_log2 <- df_fragments_log2(CompiledFragments_0_6_noFAIMS, meta_0_6_noFAIMS)

meta_0_6_FAIMS <- meta_replicates(CompiledFragments_0_6_FAIMS)
df_0_6_FAIMS_log2 <- df_fragments_log2(CompiledFragments_0_6_FAIMS, meta_0_6_FAIMS)


meta_0_8_noFAIMS <- meta_replicates(CompiledFragments_0_8_noFAIMS)
df_0_8_noFAIMS_log2 <- df_fragments_log2(CompiledFragments_0_8_noFAIMS, meta_0_8_noFAIMS)

meta_0_8_FAIMS <- meta_replicates(CompiledFragments_0_8_FAIMS)
df_0_8_FAIMS_log2 <- df_fragments_log2(CompiledFragments_0_8_FAIMS, meta_0_8_FAIMS)


meta_1_noFAIMS <- meta_replicates(CompiledFragments_1_noFAIMS)
df_1_noFAIMS_log2 <- df_fragments_log2(CompiledFragments_1_noFAIMS,meta_1_noFAIMS)

meta_1_FAIMS <- meta_replicates(CompiledFragments_1_FAIMS)
df_1_FAIMS_log2 <- df_fragments_log2(CompiledFragments_1_FAIMS, meta_1_FAIMS)


meta_1_3_noFAIMS <- meta_replicates(CompiledFragments_1_3_noFAIMS)
df_1_3_noFAIMS_log2 <- df_fragments_log2(CompiledFragments_1_3_noFAIMS, meta_1_3_noFAIMS)

meta_1_3_FAIMS <- meta_replicates(CompiledFragments_1_3_FAIMS)
df_1_3_FAIMS_log2 <- df_fragments_log2(CompiledFragments_1_3_FAIMS, meta_1_3_FAIMS)


meta_2_noFAIMS <- meta_replicates(CompiledFragments_2_noFAIMS)
df_2_noFAIMS_log2 <- df_fragments_log2(CompiledFragments_2_noFAIMS, meta_2_noFAIMS)

meta_2_FAIMS <- meta_replicates(CompiledFragments_2_FAIMS)
df_2_FAIMS_log2 <- df_fragments_log2(CompiledFragments_2_FAIMS, meta_2_FAIMS)

#### Annotations for pheatmaps ####

#x = df_0_4_noFAIMS_log2
#m = meta_0_4_noFAIMS
replicate_annotation <- function(x,m){
  annotation_res = data.frame(Resolution = m$Resolution)
  rownames(annotation_res) <- colnames(x)
  
  annotation_MI = data.frame(MaxInjection = m$MaxInjection)
  rownames(annotation_MI) <- colnames(x)
  
  annotation = cbind(annotation_res, annotation_MI)
  
  return(annotation)
}


annotation_0_4_noFAIMS <- replicate_annotation(df_0_4_noFAIMS_log2, meta_0_4_noFAIMS)
annotation_0_6_noFAIMS <- replicate_annotation(df_0_6_noFAIMS_log2, meta_0_6_noFAIMS)
annotation_0_8_noFAIMS <- replicate_annotation(df_0_4_noFAIMS_log2, meta_0_8_noFAIMS)
annotation_1_noFAIMS <- replicate_annotation(df_1_noFAIMS_log2, meta_1_noFAIMS)
annotation_1_3_noFAIMS <- replicate_annotation(df_1_3_noFAIMS_log2, meta_1_3_noFAIMS)
annotation_2_noFAIMS <- replicate_annotation(df_2_noFAIMS_log2, meta_2_noFAIMS)

annotation_0_4_FAIMS <- replicate_annotation(df_0_4_FAIMS_log2, meta_0_4_FAIMS)
annotation_0_6_FAIMS <- replicate_annotation(df_0_6_FAIMS_log2, meta_0_6_FAIMS)
annotation_0_8_FAIMS <- replicate_annotation(df_0_8_FAIMS_log2, meta_0_8_FAIMS)
annotation_1_FAIMS <- replicate_annotation(df_1_FAIMS_log2, meta_1_FAIMS)
annotation_1_3_FAIMS <- replicate_annotation(df_1_3_FAIMS_log2, meta_1_3_FAIMS)
annotation_2_FAIMS <- replicate_annotation(df_2_FAIMS_log2, meta_2_FAIMS)

# my_colour = list(
#   MaxInjection = c(`0.502` = "#FFF19C", `1` = "#FEC879", `5` = "#F8976B" ),
#   Resolution = c(`120K` = "#AEAEAE", `240K` = "#7C91AB", `500K` = "#716F9A")
#   #Resolution = c(`120K` = "#E2C044", `240K` = "#DB7C26", `500K` = "#DB5A42")
# )




Experiment_meta <- data.frame(IsolationWidth = c("0_4","0_6","0_8","1","1_3","2")) 
IW_sample_noFAIMS = vector()
  for(i in 1:length(Experiment_meta$IsolationWidth)){
    IW_sample_noFAIMS = append(IW_sample_noFAIMS, paste0("df_",Experiment_meta$IsolationWidth[i],"_noFAIMS_log2"))
  }
IW_sample_FAIMS = vector()
for(i in 1:length(Experiment_meta$IsolationWidth)){
  IW_sample_FAIMS = append(IW_sample_FAIMS, paste0("df_",Experiment_meta$IsolationWidth[i],"_FAIMS_log2"))
}
Experiment_meta$FAIMS <- IW_sample_FAIMS

noFAIMS.list <- list(df_0_4_noFAIMS,df_0_6_noFAIMS_log2,df_0_8_noFAIMS_log2,df_1_noFAIMS_log2, df_1_3_noFAIMS_log2, df_2_noFAIMS_log2)
FAIMS.list <- list(df_0_4_FAIMS_log2,df_0_6_FAIMS_log2,df_0_8_FAIMS_log2, df_1_FAIMS_log2, df_1_3_FAIMS_log2, df_2_FAIMS_log2)


my_colour = list(
  #MaxInjection = c(`0.502` = "#CEE3DD", `1` = "#A7C2A4", `5` = "#73A08A"),
  MaxInjection = c(`0.502` = "#F7F7F7", `1` = "#D0D3CA", `5` = "#A3A79D"),
  Resolution = c(`120K` = "#FFEC77", `240K` = "#FFBE5E", `500K` = "#DB664F")
  #Resolution = c(`120K` = "#E2C044", `240K` = "#DB7C26", `500K` = "#DB5A42")
)

cols <- brewer.pal(3, "YlGnBu")
scaleRYG <- colorRampPalette(cols)(20)
#scaleRYG <- colorRampPalette(c("#F5F5F5","#E6F598","#78A073"))(20)

for(i in 1:length(noFAIMS.list)){
  
  pheatmap(noFAIMS.list[[i]],
           #annotation = annotation_res,
           annotation_colors = my_colour,
           annotation_col = annotation,
           color = scaleRYG,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           gaps_col =c(9,18,27),
           gaps_row = c(1),
           show_colnames = F,
           main = "Isolation Width of 0.4 Da - Parameter Comparisons Infusions")
}


pheatmap(df_0_4_noFAIMS_log2,
         #annotation = annotation_res,
         annotation_colors = my_colour,
         annotation_col = annotation_0_4_noFAIMS,
         color = scaleRYG,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         gaps_col =c(9,18,27),
         gaps_row = c(1),
         show_colnames = F,
         main = "Isolation Width of 0.4 Da - Parameter Comparisons Infusions")

pheatmap(df_0_4_FAIMS_log2,
         #annotation = annotation_res,
         annotation_colors = my_colour,
         annotation_col = annotation_0_4_FAIMS,
         color = scaleRYG,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         gaps_col =c(9,18,27),
         gaps_row = c(1),
         show_colnames = F,
         main = "Isolation Width of 0.4 Da with FAIMS - Parameter Comparisons Infusions") 

pheatmap(df_0_6_noFAIMS_log2,
         #annotation = annotation_res,
         annotation_colors = my_colour,
         annotation_col = annotation_0_6_noFAIMS,
         color = scaleRYG,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         gaps_col =c(9,18,27),
         gaps_row = c(1),
         show_colnames = F,
         main = "Isolation Width of 0.6 Da - Parameter Comparisons Infusions")

pheatmap(df_0_6_FAIMS_log2,
         #annotation = annotation_res,
         annotation_colors = my_colour,
         annotation_col = annotation_0_6_FAIMS,
         color = scaleRYG,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         gaps_col =c(9,18,27),
         gaps_row = c(1),
         show_colnames = F,
         main = "Isolation Width of 0.6 Da with FAIMS - Parameter Comparisons Infusions") 


pheatmap(df_0_8_noFAIMS_log2,
         #annotation = annotation_res,
         annotation_colors = my_colour,
         annotation_col = annotation_0_8_noFAIMS,
         color = scaleRYG,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         gaps_col =c(9,18,27),
         gaps_row = c(1),
         show_colnames = F,
         main = "Isolation Width of 0.8 Da - Parameter Comparisons Infusions")

pheatmap(df_0_8_FAIMS_log2,
         #annotation = annotation_res,
         annotation_colors = my_colour,
         annotation_col = annotation_0_8_FAIMS,
         color = scaleRYG,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         gaps_col =c(9,18,27),
         gaps_row = c(1),
         show_colnames = F,
         main = "Isolation Width of 0.8 Da with FAIMS - Parameter Comparisons Infusions") 


pheatmap(df_1_noFAIMS_log2,
         #annotation = annotation_res,
         annotation_colors = my_colour,
         annotation_col = annotation_1_noFAIMS,
         color = scaleRYG,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         gaps_col =c(9,18,27),
         gaps_row = c(1),
         show_colnames = F,
         main = "Isolation Width of 1 Da - Parameter Comparisons Infusions")

pheatmap(df_1_FAIMS_log2,
         #annotation = annotation_res,
         annotation_colors = my_colour,
         annotation_col = annotation_1_FAIMS,
         color = scaleRYG,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         gaps_col =c(9,18,27),
         gaps_row = c(1),
         show_colnames = F,
         main = "Isolation Width of 1 Da with FAIMS - Parameter Comparisons Infusions") 

pheatmap(df_1_3_noFAIMS_log2,
         #annotation = annotation_res,
         annotation_colors = my_colour,
         annotation_col = annotation_1_3_noFAIMS,
         color = scaleRYG,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         gaps_col =c(9,18,27),
         gaps_row = c(1),
         show_colnames = F,
         main = "Isolation Width of 1.3 Da - Parameter Comparisons Infusions")

pheatmap(df_1_3_FAIMS_log2,
         #annotation = annotation_res,
         annotation_colors = my_colour,
         annotation_col = annotation_1_3_FAIMS,
         color = scaleRYG,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         gaps_col =c(9,18,27),
         gaps_row = c(1),
         show_colnames = F,
         main = "Isolation Width of 1.3 Da with FAIMS - Parameter Comparisons Infusions") 

pheatmap(df_2_noFAIMS_log2,
         #annotation = annotation_res,
         annotation_colors = my_colour,
         annotation_col = annotation_2_noFAIMS,
         color = scaleRYG,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         gaps_col =c(9,18,27),
         gaps_row = c(1),
         show_colnames = F,
         main = "Isolation Width of 2 Da - Parameter Comparisons Infusions")

pheatmap(df_2_FAIMS_log2,
         #annotation = annotation_res,
         annotation_colors = my_colour,
         annotation_col = annotation_2_FAIMS,
         color = scaleRYG,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         gaps_col =c(9,18,27),
         gaps_row = c(1),
         show_colnames = F,
         main = "Isolation Width of 2 Da with FAIMS - Parameter Comparisons Infusions") 

replicate_heatmap(df_0_4_noFAIMS_log2, meta_0_8_noFAIMS)

# pheatmap(df_0_4_noFAIMS_log2,
#          #annotation = annotation_res,
#          annotation_col = annotation_MI,
#          color = scaleRYG,
#          cluster_rows = FALSE,
#          cluster_cols = FALSE,
#          gaps_col =seq(from = 3, to = 24, by = 3),
#          #gaps_row = c(1,10),
#          main = "Isolation Width of 0.4 Da - Parameter Comparisons Infusions")

pheatmap(df_0_4_FAIMS_log2,
         color = scaleRYG,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         main = "Isolation Width of 0.4 Da - Parameter Comparisons Infusions With FAIMS")

pheatmap(df_0_6_noFAIMS_log2,
         color = scaleRYG,
         cluster_rows = FALSE,
         cluster_cols = FALSE)

pheatmap(df_0_6_FAIMS_log2,
         color = scaleRYG,
         cluster_rows = FALSE,
         cluster_cols = FALSE)

pheatmap(df_0_8_noFAIMS_log2,
         color = scaleRYG,
         cluster_rows = FALSE,
         cluster_cols = FALSE)

pheatmap(df_0_8_FAIMS_log2,
         color = scaleRYG,
         cluster_rows = FALSE,
         cluster_cols = FALSE)

pheatmap(df_1_noFAIMS_log2,
         color = scaleRYG,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         main = "Isolation Width of 1 Da - Parameter Comparisons Infusions")

pheatmap(df_1_FAIMS_log2,
         color = scaleRYG,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         main = "Isolation Width of 1 Da - Parameter Comparisons Infusions with FAIMS")

pheatmap(df_1_3_noFAIMS_log2,
         color = scaleRYG,
         cluster_rows = FALSE,
         cluster_cols = FALSE)

pheatmap(df_2_FAIMS_log2,
         color = scaleRYG,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         main = "Isolation Width of 2 Da - Parameter Comparisons Infusions")

pheatmap(df_2_noFAIMS_log2,
         color = scaleRYG,
         cluster_rows = FALSE,
         cluster_cols = FALSE)

# To combine noFaims and FAIMs samples into one heatmap
fragments <- rownames(df_0_4_FAIMS_log2)
fragment_FAIMS = vector()
for(i in fragments){fragment_FAIMS = append(fragment_FAIMS,paste0(i,'_FAIMS', collapse = " "))}

rownames(df_0_4_FAIMS_log2) <- fragment_FAIMS

df_0_4_Combined_Log2 <- rbind(df_0_4_noFAIMS_log2, df_0_4_FAIMS_log2)

pheatmap(df_0_4_Combined_Log2,
         #color = scaleRYG,
         cluster_rows = FALSE,
         cluster_cols = FALSE)

#Clean up heat map so that MaxInjection is ordered. 
colnames(df_0_4_noFAIMS_log2)

pdf("H:/Projects/Proteomics/Zymomona/FAIMS/Figures/FromR/AIEIVDQALDR_combined_0_4_log2_heatmap.pdf", width = 3, height = 7,useDingbats = FALSE)
pheatmap(Combined_0_4_Log2,
         annotation_colors = my_colour,
         annotation_col = unique_annotation_0_4_noFAIMS,
         #annotation_row = fragment_annotation,
         color = scaleRYG,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         gaps_col =c(3,6,9),
         gaps_row = c(1,19,20),
         show_colnames = F,
         legend = F,
         annotation_legend = F,
         main = "AIEIVDQALDR IW of 0.4 Da")
dev.off()

