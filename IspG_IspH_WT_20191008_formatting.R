library(installr) # for is.empty function
#BiocManager::install("pcaMethods")
library(pcaMethods)
library(ggplot2)
library(plotly)
library(reshape2)

#browseVignettes("pcaMethods")

#### Upload Data ####

palette(c("#5ABCB2", "#D6F8D6","#D6F8D6" ,"mediumorchid2","darkgoldenrod1", 
          "hotpink3", "red2", "sienna2","slategray4","mediumturquoise", 
          "deepskyblue", "orangered", "midnightblue"))

# Read in MaxQuant output proteingroups, peptides as well as DRB imputed dataframe in Data_Formatting_Functions
load("P:/EAT_20190926_Zymomona/ISPH_ISPG_OE_Mehmet/EclipseRuns/combined_ZM4tag/txt/OEZM4tag_45minLCMSMS_Eclipse_AbundanceImputed_20200420.RData")

# load 60 min LC-MS/MS data collected on Boudicca? confirm.
proteinGroups_Zymo_Eclipse60min <- read.csv("P:/EAT_20190926_Zymomona/ISPH_ISPG_OE_Mehmet/txt_60min_combinedSearch/proteinGroups.csv", 
                                            header = TRUE, sep = ",", stringsAsFactors = FALSE)

rownames(merged_drb_impute) <- merged_drb_impute$Majority.protein.ID

#### ISPG Protein #####
#subset the protein groups that are associated to the ISPG gene

ISPG_peptides <- peptides_Zymo_Eclipse45min_ZM4tag[grep("ZMO0180", peptides_Zymo_Eclipse45min_ZM4tag$Proteins),]

ISPG_reformat_uniprot <- function(gene,strain,df){
  ISPG_data <- df[grep(gene, df$Fasta.headers),]
  ISPG_lfq <- ISPG_data[grep(strain, ISPG_data$Fasta.headers),c(1,grep("LFQ",colnames(ISPG_data)))]
  Sample <- sub("LFQ.intensity.", "", colnames(ISPG_lfq[,-1]))
  Intensity <- unlist(ISPG_lfq[2:ncol(ISPG_lfq)])
  ISPG_lfq <- data.frame("Sample" = Sample, "Intensity" = Intensity)
  ISPG_lfq$Type <- c("ISPG","ISPG","ISPG","ISPH","ISPH","ISPH","WT","WT","WT")
}


#df <- merged_drb_impute
#df_peptide <- ISIPG_peptides
ISPG_reformat_ZM4 <- function(gene,df, peptide,df_peptide){
  ISPG_data <- df[grep(gene, df$Protein.IDs),]
  ISPG_lfq <- ISPG_data[, c(1,grep("LFQ", colnames(ISPG_data)))]

  for(i in 1:length(peptide))
  {
    ISPG_peptide_lfq <- df_peptide[which(df_peptide$Sequence == peptide[i]),grep("Intensity", colnames(df_peptide))] 
    ISPG_peptide_lfq <- ISPG_peptide_lfq[,-grep("new", colnames(ISPG_peptide_lfq))]
    assign(paste(peptide[i],"_lfq",sep = ""), vapply(ISPG_peptide_lfq[1,-1], paste, collapse = ",", character(1L)))

  }
  
  Sample <- sub("LFQ.intensity.", "", colnames(ISPG_lfq[,-1]))
  Intensity <- unlist(ISPG_lfq[2:ncol(ISPG_lfq)])
  ISPG_lfq <- data.frame("Sample" = Sample, "Intensity" = Intensity)
  ISPG_lfq$VSLSADPEQEVR <- as.numeric(VSLSADPEQEVR_lfq)
  ISPG_lfq$ETDIGVTGGGQGK <- as.numeric(ETDIGVTGGGQGK_lfq)
  ISPG_lfq$AVDCPLHLGITEAGGLIGGTVK <- as.numeric(AVDCPLHLGITEAGGLIGGTVK_lfq)
  ISPG_lfq$Type <- c("ISPG","ISPG","ISPG","ISPH","ISPH","ISPH","WT","WT","WT")
  
  v = vector()
  WT_lfq <- ISPG_lfq[grep("WT",rownames(ISPG_lfq)),2]
  
    for (j in 1:nrow(ISPG_lfq)) {
      for (i in 1:3) {
      if(length(grep(i,ISPG_lfq$Sample[j])) != 0){
        v <- append(v,ISPG_lfq$Intensity[j]/WT_lfq[i])
        print(v)
        }
      }
    }
  ISPG_lfq$FC <- v
  return(ISPG_lfq)
}

# df <- ISPG_lfq
# colnum will be 3,4,5
peptide_FC <- function(colnum,df){
  v = numeric(0)
  WT_lfq <- as.numeric(df[grep("WT",rownames(df)),colnum])
  
  for (j in 1:nrow(df)) {
    for (i in 1:3) {
      if(length(grep(i,df$Sample[j])) != 0){
        v <- append(v,(as.numeric(df[j,colnum])/WT_lfq[i]))
        print(v)
      }
    }
  }
  
  return(v)
}



peptide <- c("VSLSADPEQEVR", "ETDIGVTGGGQGK", "AVDCPLHLGITEAGGLIGGTVK")
# Apply function to ISPG gene "ZMO0180"
ISPG_lfq <- ISPG_reformat_ZM4("ZMO0180", merged_drb_impute, peptide, ISPG_peptides)

ISPG_lfq$VSLSADPEQEVR_FC <- peptide_FC(3,ISPG_lfq)
ISPG_lfq$ETDIGVTGGGQGK_FC <- peptide_FC(4,ISPG_lfq)
ISPG_lfq$AVDCPLHLGITEAGGLIGGTVK_FC <- peptide_FC(5,ISPG_lfq)

ISPG_lfq$VSLSADPEQEVR_FC_log2 <- log2(ISPG_lfq$VSLSADPEQEVR_FC)
ISPG_lfq$ETDIGVTGGGQGK_FC__log2 <- log2(ISPG_lfq$ETDIGVTGGGQGK_FC)
ISPG_lfq$AVDCPLHLGITEAGGLIGGTVK_FC__log2 <- log2(ISPG_lfq$AVDCPLHLGITEAGGLIGGTVK_FC)

#### ISPG Plotting ####

type.colors = as.numeric(factor(ISPG_lfq$Type))
q <- ggplot(ISPG_lfq, aes(x=Sample, y=log2(Intensity), fill=Type)) +
  geom_bar(stat = "identity")+
  scale_fill_manual(values = c("#5ABCB2", "#C3DE79","#A394DD"))+
  theme_light()
q + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(title = "Log2(ISPG protein)", 
                                                                    subtitle = "Across sample types in triplicate",
                                                                    y = "Log2(Intensity)")

type.colors = as.numeric(factor(ISPG_lfq$Type))
q <- ggplot(ISPG_lfq, aes(x=Type, y=FC, fill=Type)) +
  geom_boxplot()+
  scale_fill_manual(values = c("#5ABCB2", "#C3DE79","#A394DD"))+
  theme_light()
q + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(title = "Fold Change ISPG protein", 
                                                                    subtitle = "Over Expression (OE) Strain",
                                                                    y = "Fold Change (OE/WT)")
type.colors = as.numeric(factor(ISPG_lfq$Type))
q <- ggplot(ISPG_lfq, aes(x=Type, y=log2(Intensity), fill=Type)) +
  geom_boxplot()+
  geom_point(data = ISPG_lfq, aes(x= Type, y= log2(Intensity)), size = 3, color = "#6D696F", alpha = 0.5) +
  scale_fill_manual(values = c("#5ABCB2", "#C3DE79","#A394DD"))+
  theme_light()
q + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(title = "Log2(Abundance) ISPG protein", 
                                                                    subtitle = "Across sample types in triplicate",
                                                                    y = "Log2(Intensity)")

#### ISPG Peptides ####
#subset zymo peptides for protein Q5NR50 or ZMO0180

ISPG_peptides_averages <- ISPG_peptides[,c(42,grep("IspG", colnames(ISPG_peptides)))]
sub_ISPG_peptides_avg <- ISPG_peptides_averages[,c(1,11:13)]
colnames(sub_ISPG_peptides_avg) <- c("PEP","ISPG-1", "ISPG-2", "ISPG-3")
rownames(sub_ISPG_peptides_avg) <- ISPG_peptides$Sequence

sub_ISPG_peptides_avg$Avg <- rowMeans(sub_ISPG_peptides_avg[,2:ncol(sub_ISPG_peptides_avg)])
sub_ISPG_peptides_avg$Log2Avg <- log2(sub_ISPG_peptides_avg$Avg)
sub_ISPG_peptides_avg$Std <- apply(sub_ISPG_peptides_avg[,2:4],1,function(x) sd(log2(x)))


v <- vector()
#subset peptides with abundance greater than log2
for(i in 1:nrow(sub_ISPG_peptides_avg)){
  if(0 %in% sub_ISPG_peptides_avg[i,2:4])
    v <- append(v,i) 
}

ISPG_peptides_abundance <- sub_ISPG_peptides_avg[-v,]
ISPG_peptides_abundance <- ISPG_peptides_abundance[order(ISPG_peptides_abundance$Log2Avg),]
ISPG_peptides_abundance$Peptide <- rownames(ISPG_peptides_abundance)

ISPG_peptides_abundance<- ISPG_peptides_abundance %>%
  mutate(Peptide = factor(Peptide, level = rownames(ISPG_peptides_abundance)))

str(ISPG_peptides_abundance)
ISPG_peptides_abundance$Peptide


long_ISPG_Lfq <- melt(ISPG_lfq, id.vars = c("Sample"))

# ggplot(ISPG_peptides_abundance, aes(x=Peptide, y= Log2Avg, )) + 
#   geom_point(aes(size = PEP)) +
#   scale_color_continuous(low="blue", high = "red", breaks=c(2.59e-238, 0.090302)) +
#   scale_size(breaks = c(2.59e-238,0.090302), range= c(1,15)) +
#   theme_light()+
#   labs(title = "Average Peptide Abundances" ,subtitle = "IspG Over Expressed Samples", y = "Log2(Intensity)", x ="Peptides") +
#   coord_flip()
# #   theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(ISPG_peptides_abundance, aes(x=Peptide, y= Log2Avg, )) +
  geom_point() +
  theme_light()+
  labs(title = "Average Peptide Abundances" ,subtitle = "IspG Over Expressed Samples", y = "Log2(Intensity)", x ="Peptides") +
  coord_flip()
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
 
ggplot(ISPG_peptides_abundance, aes(x=Peptide, y= Log2Avg)) + 
  geom_bar(stat = "identity", color = "#6D696F" ,  fill = "#7C9E3D", alpha = 0.6, position = position_dodge())+
  geom_errorbar(aes(ymin=Log2Avg, ymax=Log2Avg+Std) , width = .2, color= "#6D696F", position = position_dodge(0.052)) +
  geom_point(data = ISPG_peptides_abundance, aes(x=Peptide, y = log2(ISPG_peptides_abundance$`ISPG-1`)), size = 3, color = "#6D696F", alpha = 0.5) +
  geom_point(data = ISPG_peptides_abundance, aes(x= Peptide, y= log2(ISPG_peptides_abundance$`ISPG-2`)), size = 3, color = "#6D696F", alpha = 0.5) +
  geom_point(data = ISPG_peptides_abundance, aes(x= Peptide, y= log2(ISPG_peptides_abundance$`ISPG-3`)), size = 3, color = "#6D696F", alpha = 0.5) +
  theme_light()+
  labs(title = "Peptide Abundances" ,subtitle = "IspG Over Expressed Samples", y = "Log2(Intensity)", x ="Peptides") +
  coord_flip()
# theme(axis.text.x = element_text(angle = 45, hjust = 1)) 



type.colors = as.numeric(factor(ISPG_lfq$Type))
q <- ggplot(ISPG_lfq, aes(x=Type, y=VSLSADPEQEVR_FC, fill=Type)) +
  geom_boxplot()+
  geom_boxplot(data = ISPG_lfq, aes(x=Type, y = ETDIGVTGGGQGK_FC, fill = Type))+
  scale_fill_manual(values = c("#5ABCB2", "#C3DE79","#A394DD"))+
  theme_light()
q + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(title = "Fold Change ISPG protein", 
                                                                    subtitle = "Over Expression (OE) Strain",
                                                                    y = "Fold Change (OE/WT)")
#### ISPH Protein #####

#gene <- "GN=ispH"
#strain <- "ZM4"
#df <- merged_drb_impute
ISPH_reformat_uniprot <- function(gene,strain,df){
  ISPH_data <- df[grep(gene, df$Fasta.headers),]
  ISPH_lfq <- ISPH_data[grep(strain, ISPH_data$Fasta.headers),c(1,grep("LFQ",colnames(ISPH_data)))]
  Sample <- sub("LFQ.intensity.", "", colnames(ISPH_lfq[,-1]))
  Intensity <- unlist(ISPH_lfq[2:ncol(ISPH_lfq)])
  ISPH_lfq <- data.frame("Sample" = Sample, "Intensity" = Intensity)
  ISPH_lfq$Type <- c("ISPG","ISPG","ISPG","ISPH","ISPH","ISPH","WT","WT","WT")
}


#df <- merged_drb_impute
ISPH_reformat_ZM4 <- function(gene,df){
  ISPH_data <- df[grep(gene, df$Protein.IDs),]
  ISPH_lfq <- ISPH_data[, c(1,grep("LFQ", colnames(ISPH_data)))]
  Sample <- sub("LFQ.intensity.", "", colnames(ISPH_lfq[,-1]))
  Intensity <- unlist(ISPH_lfq[2:ncol(ISPH_lfq)])
  ISPH_lfq <- data.frame("Sample" = Sample, "Intensity" = Intensity)
  ISPH_lfq$Type <- c("ISPG","ISPG","ISPG","ISPH","ISPH","ISPH","WT","WT","WT")
  
  v = vector()
  WT_lfq <- ISPG_lfq[grep("WT",rownames(ISPH_lfq)),2]
  
  for (j in 1:nrow(ISPH_lfq)) {
    for (i in 1:3) {
      if(length(grep(i,ISPH_lfq$Sample[j])) != 0){
        v <- append(v,ISPH_lfq$Intensity[j]/WT_lfq[i])
        print(v)
      }
    }
  }
  ISPH_lfq$FC <- v
  return(ISPH_lfq)
}

# Apply function to ISPG gene "ZMO0875"
ISPH_lfq <- ISPH_reformat_ZM4("ZMO0180", merged_drb_impute)

type.colors = as.numeric(factor(ISPH_lfq$Type))
q <- ggplot(ISPH_lfq, aes(x=Sample, y=Intensity, fill=Type)) +
  geom_bar(stat = "identity")+
  scale_fill_manual(values = c("#5ABCB2", "#C3DE79","#A394DD"))+
  theme_light()
q + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(title = "Log2(ISPH protein)", 
                                                                    subtitle = "Across sample types in triplicate",
                                                                    y = "Log2(Intensity)")

type.colors = as.numeric(factor(ISPH_lfq$Type))
q <- ggplot(ISPH_lfq, aes(x=Type, y=FC, fill=Type)) +
  geom_boxplot()+
  scale_fill_manual(values = c("#5ABCB2", "#C3DE79","#A394DD"))+
  theme_light()
q + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(title = "Fold Change ISPH protein", 
                                                                    subtitle = "Over Expression (OE) Strain",
                                                                    y = "Fold Change (OE/WT)")
type.colors = as.numeric(factor(ISPH_lfq$Type))
q <- ggplot(ISPH_lfq, aes(x=Type, y=log2(Intensity), fill=Type)) +
  geom_boxplot()+
  geom_point(data = ISPH_lfq, aes(x= Type, y= log2(Intensity)), size = 3, color = "#6D696F", alpha = 0.5) +
  scale_fill_manual(values = c("#5ABCB2", "#C3DE79","#A394DD"))+
  theme_light()
q + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(title = "Log2(Abundance) ISPH protein", 
                                                                    subtitle = "Across sample types in triplicate",
                                                                    y = "Log2(Intensity)")
#### ISPH Peptide ####
#subset zymo peptides for protein Q5NR50 or "ZMO0875"
ISPH_peptides <- peptides_Zymo_Eclipse45min_ZM4tag[grep("ZMO0875", peptides_Zymo_Eclipse45min_ZM4tag$Proteins),]

ISPH_peptides_averages <- ISPH_peptides[,c(42,grep("IspH", colnames(ISPH_peptides)))]
sub_ISPH_peptides_avg <- ISPH_peptides_averages[,c(1,8:10)]
colnames(sub_ISPH_peptides_avg) <- c("PEP","ISPH-1", "ISPH-2", "ISPH-3")
rownames(sub_ISPH_peptides_avg) <- ISPH_peptides$Sequence

sub_ISPH_peptides_avg$Avg <- rowMeans(sub_ISPH_peptides_avg[,2:ncol(sub_ISPH_peptides_avg)])
sub_ISPH_peptides_avg$Log2Avg <- log2(sub_ISPH_peptides_avg$Avg)
sub_ISPH_peptides_avg$Std <- apply(sub_ISPH_peptides_avg[,2:4],1,function(x) sd(log2(x)))


v <- vector()
#subset peptides with abundance greater than log2
for(i in 1:nrow(sub_ISPH_peptides_avg)){
  if(0 %in% sub_ISPH_peptides_avg[i,2:4])
    v <- append(v,i) 

}

ISPH_peptides_abundance <- sub_ISPH_peptides_avg[-v,]
ISPH_peptides_abundance <- ISPH_peptides_abundance[order(ISPH_peptides_abundance$Log2Avg),]
ISPH_peptides_abundance$Peptide <- rownames(ISPH_peptides_abundance)

ISPH_peptides_abundance<- ISPH_peptides_abundance %>%
  mutate(Peptide = factor(Peptide, level = rownames(ISPH_peptides_abundance)))

str(ISPH_peptides_abundance)
ISPG_peptides_abundance$Peptide

ggplot(ISPH_peptides_abundance, aes(x=Peptide, y= Log2Avg)) + 
  geom_point() +
  theme_light()+
  labs(title = "Average Peptide Abundances" ,subtitle = "IspH Over Expressed Samples", y = "Log2(Intensity)", x ="Peptides") +
  coord_flip()
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggplot(ISPH_peptides_abundance, aes(x=Peptide, y= Log2Avg)) + 
  geom_bar(stat = "identity", color = "#6D696F" ,  fill = "#E0D68A", alpha = 0.6, position = position_dodge())+
  geom_errorbar(aes(ymin=Log2Avg, ymax=Log2Avg+Std) , width = .5, color = "#6D696F", position = position_dodge(0.052)) +
  geom_point(data = ISPH_peptides_abundance, aes(x=Peptide, y = log2(ISPH_peptides_abundance$`ISPH-1`)), color = "#6D696F",  size = 4, alpha = 0.5) +
  geom_point(data = ISPH_peptides_abundance, aes(x= Peptide, y= log2(ISPH_peptides_abundance$`ISPH-2`)), color = "#6D696F", size = 4, alpha = 0.5) +
  geom_point(data = ISPH_peptides_abundance, aes(x= Peptide, y= log2(ISPH_peptides_abundance$`ISPH-3`)), color = "#6D696F", size = 4, alpha = 0.5) +
  theme_light()+
  labs(title = "Peptide Abundances" ,subtitle = "IspH Over Expressed Samples", y = "Log2(Intensity)", x ="Peptides") +
  coord_flip()
# theme(axis.text.x = element_text(angle = 45, hjust = 1)) 


#### DXS2 peptides #####
#subset the protein groups that are associated to the ISPH gene
DXS2_data <- proteinGroups_Zymo_Eclipse60min[grep("GN=DXS2", proteinGroups_Zymo_Eclipse60min$Fasta.headers),]

ISPH_lfq <- ISPH_data[,c(1,grep("LFQ",colnames(ISPH_data)))]

#### ZM4 proteins ####
proteinGroups_ZM4 <- merged_drb_impute[grep("ZM4",merged_drb_impute$Fasta.headers),]
proteinGroups_ZM4 <- proteinGroups_ZM4[order(-proteinGroups_ZM4$Average.IspG),]
proteinGroups_ZM4 <- proteinGroups_ZM4[order(-proteinGroups_ZM4$Average.IspH),]
proteinGroups_ZM4 <- proteinGroups_ZM4[order(-proteinGroups_ZM4$Average.WT),]

grep("ISPG", proteinGroups_ZM4$Fasta.headers)
proteinGroups_ZM4[162,] #ISPH

p <- ggplot(merged_drb_impute, aes(y=Average.IspH, x = Average.IspG)) +
  geom_point()

ggplotly(p)

merged_drb_impute[grep(32.16088, merged_drb_impute$Average.IspG),]

p <- ggplot(merged_drb_impute, aes(y=Average.IspH, x = Average.WT)) +
  geom_point()

ggplotly(p)


p <- ggplot(merged_drb_impute, aes(y=Average.IspG, x = Average.WT)) +
  geom_point()

ggplotly(p)

merged_drb_impute[grep(37.704, merged_drb_impute$Average.IspG),]

ggplot(merged_drb_impute, aes(x=pvalues)) + 
  geom_histogram(binwidth = .05, color="black", fill="white")




df <- data.frame(matrix(unlist(drb_merged[,-1]), nrow=9, byrow=T),stringsAsFactors=FALSE)
df <- t(df)
colnames(df) <- colnames(drb_merged[,-1])
rownames(df) <- drb_merged[,1]

pvals=apply(df,1,function(x) {t.test(x[1:3],x[7:9])$p.value})

merged_drb_impute$pvalues <- pvals

var(df[1:3])
var(df[7:9])

hist(merged_drb_impute[16:24], breaks = 50)
boxplot(log2(clean_df[,1:ncol(clean_df)]))

matrix_clean_df <- matrix(unlist(clean_df),ncol = 9, byrow = FALSE)
colnames(matrix_clean_df) <- colnames(clean_df)
rownames(matrix_clean_df) <- rownames(clean_df)

# no missing values allowed #resPCA <- pca(matrix_clean_df, method="svd", center=FALSE, nPcs=5)
resPPCA <- pca(matrix_clean_df, method="ppca", center=FALSE, nPcs=5)
resBPCA <- pca(matrix_clean_df, method="bpca", center=FALSE, nPcs=5)
resSVDI <- pca(matrix_clean_df, method="svdImpute", center=FALSE, nPcs=5)
resNipals <- pca(matrix_clean_df, method="nipals", center=FALSE, nPcs=5)
resNLPCA <- pca(matrix_clean_df, method="nlpca", center=FALSE, nPcs=5, maxSteps=300)
