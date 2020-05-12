library(pcaMethods)
library(ggplot2)
library(pheatmap)
library(plotly)
library(reshape2)
library(RColorBrewer)
library(plyr)

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
df_fragments_log2 <- function(x){
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
  
  return(df_log2)
  
}



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
  
  meta[,2] <- sub('MI', '', meta[,2])
  meta[,2] <- sub(502, 0.502, meta[,2])
  meta[,2] <- as.numeric(meta[,2])

}

df_0_4_noFAIMS_log2 <- df_fragments_log2(CompiledFragments_0_4_noFAIMS)
df_0_4_FAIMS_log2 <- df_fragments_log2(CompiledFragments_0_4_FAIMS)

df_0_6_noFAIMS_log2 <- df_fragments_log2(CompiledFragments_0_6_noFAIMS)
df_0_6_FAIMS_log2 <- df_fragments_log2(CompiledFragments_0_6_FAIMS)

df_0_8_noFAIMS_log2 <- df_fragments_log2(CompiledFragments_0_8_noFAIMS)
df_0_8_FAIMS_log2 <- df_fragments_log2(CompiledFragments_0_8_FAIMS)

df_1_noFAIMS_log2 <- df_fragments_log2(CompiledFragments_1_noFAIMS)
df_1_FAIMS_log2 <- df_fragments_log2(CompiledFragments_1_FAIMS)

df_1_3_noFAIMS_log2 <- df_fragments_log2(CompiledFragments_1_3_noFAIMS)
df_1_3_FAIMS_log2 <- df_fragments_log2(CompiledFragments_1_3_FAIMS)

df_2_noFAIMS_log2 <- df_fragments_log2(CompiledFragments_2_noFAIMS)
df_2_FAIMS_log2 <- df_fragments_log2(CompiledFragments_2_FAIMS)

cols <- brewer.pal(3, "YlGnBu")
scaleRYG <- colorRampPalette(cols)(20)

pheatmap(df_0_4_noFAIMS_log2,
         color = scaleRYG,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         main = "Isolation Width of 0.4 Da - Parameter Comparisons Infusions")

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
