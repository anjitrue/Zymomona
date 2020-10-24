library(installr) # for is.empty function

#### Format Data Functions ####

#### Upload Data ####

palette(c("#5ABCB2", "#D6F8D6","#D6F8D6" ,"mediumorchid2","darkgoldenrod1", 
          "hotpink3", "red2", "sienna2","slategray4","mediumturquoise", 
          "deepskyblue", "orangered", "midnightblue"))

# Read in MaxQuant output proteingroups, I filtered out the contaminants and the fasta header overspill

# data.frame of Eclipse 45 min runs using Eclipse02. Digestion was performed in Septemeber -  induction for 180 min and the current
# data is from time point 120 for WT, IspH, and IspG.

# Note I manually removed Contaminants - Only identified by site, Reverse, Potential Contaminant
proteinGroups_Zymo_Eclipse45min <- read.csv("E:/EAT/Projects/Proteomics/Infusion_FAIMS/txt_DDA_45minHighResMS2_120TimePoint/proteinGroups.csv", 
                                            header = TRUE, sep = ",", stringsAsFactors = FALSE)

peptides_Zymo_Eclipse45min <- read.csv("E:/EAT/Projects/Proteomics/Infusion_FAIMS/txt_DDA_45minHighResMS2_120TimePoint/peptides.csv", 
                                       header = TRUE, sep = ",",stringsAsFactors = FALSE) 

# Potentially impute later?
# DRB_imputed_data <- read.csv("P:/EAT_20190926_Zymomona/ISPH_ISPG_OE_Mehmet/EclipseRuns/combined_ZM4tag/txt/Zymo_OE_ZM4tag_proteinGroups.csv")


#### Data Formating #####

# subset specific columns and bind function
subsetLFQ <- function(x){
  y<- x[,which(names(x) %in% c("Protein.IDs", 
                               "Majority.protein.IDs",
                               "Fasta.headers",
                               "Peptides",
                               "Razor...unique.peptides",
                               "Unique.peptides",
                               "Sequence.coverage....",
                               "Unique.sequence.coverage....",
                               "Mol..weight..kDa.",
                               "Q.value",
                               "Score",
                               "Intensity",
                               "MS.MS.count",
                               "id"))]
  z <- x[,grep("LFQ.intensity.",names(x))]
  z[z == 0] <- NA
  x <- cbind(y,z)
  return(x)
}

# Perform subset function
subset_df <- subsetLFQ(proteinGroups_Zymo_Eclipse45min)


# convert Protein.IDs to characters 
subset_df$Protein.IDs <- as.character(subset_df$Protein.IDs)
# convert Majority.protein.IDs to characters
subset_df$Majority.protein.IDs <- as.character(subset_df$Majority.protein.IDs)

# use id as identifier in row names
rownames(subset_df) <- subset_df$id

# Isolate only the first protein in protein group for Majority.protein.IDs and Protein.IDs
majorityProtein <- lapply(strsplit(subset_df$Majority.protein.IDs, ";"), '[', 1)
protein.ID <- lapply(strsplit(subset_df$Protein.IDs, ";"), '[', 1)

# Replace with reduced array
subset_df$Protein.IDs <- protein.ID
subset_df$Majority.protein.IDs <- majorityProtein


# remove protein groups that are missing over 50% of the data across samples

# *** corrected function below
# remove.features.50percentcuttoff <- function (x) {
#   
#   # Calculate how many missing values per feature
#   features.missing = rowMeans(is.na(x)) 
#   print(paste0("Number of protein groups that have over 50% missing measurements: ",sum(features.missing > 0.50))) 
#   features.missing.50more = rownames(x)[features.missing > 0.50] 
#   
#   # create a vector of the features that meet criteria
#   keep.features = which(features.missing <= 0.50) 
#   print(paste0("Protein groups that pass the 50% filteration: ", length(keep.features)))
#   names(keep.features) = keep.features 
#   
#   # create a vector of the features that will be removed from data.frame
#   remove.features = which(features.missing > 0.50)
#   print(paste0("Number of protein groups removed from dataset: ", length(remove.features)))
#   names(remove.features) = remove.features
#   
#   # If there isn't any features to remove keep all features
#   if(is.empty(remove.features) == TRUE ){
#     filtered = x[which(rownames(x) %in% keep.features),]
#     # otherwise filter out the features that contain over 50% of data missing
#   } else{
#     filtered = x[-which(rownames(x) %in% remove.features),]
#   }
#   return(filtered)
# }

remove.features.50percentcuttoff <- function (x) {
  
  # Calculate how many missing values per feature
  features.missing = rowMeans(is.na(x)) 
  print(paste0("Number of protein groups that have over 50% missing measurements: ",sum(features.missing > 0.50))) 
  features.missing.50more = rownames(x)[features.missing > 0.50] 
  
  # create a vector of the features that meet criteria
  keep.features = rownames(x)[features.missing <= 0.50]
  print(paste0("Protein groups that pass the 50% filteration: ", length(keep.features)))
  #names(keep.features) = keep.features 
  
  # create a vector of the features that will be removed from data.frame
  remove.features = rownames(x)[features.missing > 0.50]
  print(paste0("Number of protein groups removed from dataset: ", length(remove.features)))
  #names(remove.features) = remove.features
  
  # If there isn't any features to remove keep all features
  if(sum(features.missing) == 0 ){
    filtered = x[which(rownames(x) %in% keep.features),]
    # otherwise filter out the features that contain over 50% of data missing
  } else{
    filtered = x[-which(rownames(x) %in% remove.features),]
  }
  return(filtered)
}


# Perform feature removal on only the LFQ data only
clean_df <- remove.features.50percentcuttoff(subset_df[15:ncol(subset_df)])
# Remove Intensity from column names
colnames(clean_df) <- sub(".*intensity.", "", colnames(clean_df))

# Create a meta object with only the meta data from MaxQuant
meta <- subset_df[,1:14]

# Filtered meta data frame with the correct protein groups after 50% data cut off
filtered_meta = meta[which(rownames(meta) %in% rownames(clean_df)),]


# Bind together meta data and data
ready_to_impute <- cbind(filtered_meta,clean_df)
ready_to_impute$Protein.IDs <- vapply(ready_to_impute$Protein.IDs, paste, collapse = ", ", character(1L))
ready_to_impute$Majority.protein.IDs <- vapply(ready_to_impute$Majority.protein.IDs, paste, collapse = ", ", character(1L))
#write.csv(ready_to_impute, "P:/EAT_20190926_Zymomona/Zymo_TimeCourse/LC_RawFiles/Oct2020_Technical_Eclispe/txt_DDA_45minHighResMS2_120TimePoint/DDA_MaxQuant_proteinGroups_readyforImputation.csv")

####### import DRBimputed dataset ######
DRB_imputed_data <- read.csv("P:/EAT_20190926_Zymomona/Zymo_TimeCourse/LC_RawFiles/Oct2020_Technical_Eclispe/txt_DDA_45minHighResMS2_120TimePoint/DRBimputed_DDA_MaxQuant_proteinGroups_readyforImputation.csv")
rownames(DRB_imputed_data) <- DRB_imputed_data$Majority.protein.ID


# DRB imputed data will have meta data appended
drb_merged <- DRB_imputed_data

merged_drb_impute <- cbind(filtered_meta, drb_merged[,-1])

#### Statistics ####

#### IspG Stats ####
IspG_average <- drb_merged[,grep("ISPG", colnames(drb_merged))]
IspG_average$Average.IspG <- rowMeans(IspG_average)
IspG_average$Std.IspG <- apply(IspG_average[1:3],1,sd)


#### IspH Stats ####
IspH_average <- drb_merged[,grep("ISPH", colnames(drb_merged))]
IspH_average$Average.IspH <- rowMeans(IspH_average)
IspH_average$Std.IspH <- apply(IspH_average[1:3],1,sd)

#### WT Stats ####
WT_average <- drb_merged[,grep("WT", colnames(drb_merged))]
WT_average$Average.WT <- rowMeans(WT_average)
WT_average$Std.WT <- apply(WT_average[1:3],1,sd)

#### Fold Change ####
IspG_average$FC_IspG_WT <- IspG_average$Average.IspG-WT_average$Average.WT

IspH_average$FC_IspH_WT <- IspH_average$Average.IspH-WT_average$Average.WT


#### Append columns to imputed data set ####
merged_drb_impute$Average.IspG <- IspG_average$Average.IspG
merged_drb_impute$Average.IspH <- IspH_average$Average.IspH
merged_drb_impute$Average.WT <- WT_average$Average.WT


merged_drb_impute$FC_IspG_WT <- IspG_average$FC_IspG_WT
merged_drb_impute$FC_IspH_WT <- IspH_average$FC_IspH_WT

merged_drb_impute$Std.IspG <- IspG_average$Std.IspG
merged_drb_impute$Std.IspH <- IspH_average$Std.IspH
merged_drb_impute$Std.WT <- WT_average$Std.WT


#### t.test #####

merged_drb_impute$t.test.Ispg <- apply(cbind(IspG_average[1:3], WT_average[1:3]),1,function(x) {t.test(x[1:3], x[4:6])$p.value})
merged_drb_impute$t.test.Isph <- apply(cbind(IspH_average[1:3], WT_average[1:3]),1,function(x) {t.test(x[1:3], x[4:6])$p.value})

merged_drb_impute$FDR_corrected_Ispg <- p.adjust(merged_drb_impute$t.test.Ispg, method = "fdr")
merged_drb_impute$FDR_corrected_Isph <- p.adjust(merged_drb_impute$t.test.Isph, method = "fdr")


#### Transform lists into vectors for a complete data frame
merged_drb_impute$Protein.IDs <- vapply(merged_drb_impute$Protein.IDs, paste, collapse = ",", character(1L))
merged_drb_impute$Majority.protein.IDs <- vapply(merged_drb_impute$Majority.protein.IDs, paste, collapse = ",", character(1L))



# Write out cv
#write.csv(merged_drb_impute, "P:/EAT_20190926_Zymomona/ISPH_ISPG_OE_Mehmet/EclipseRuns/txt_Eclipse_45min/Zymo20200406_UniprotDatabase_Mean_STD_FC_proteinGroups.csv")
write.csv(ready_to_impute, "P:/EAT_20190926_Zymomona/Zymo_TimeCourse/LC_RawFiles/Oct2020_Technical_Eclispe/txt_DDA_45minHighResMS2_120TimePoint/EATsimpleStats_DRBimputed_50percentfiltered_proteinGroups_DDA_45minHighResMS2_120TimePoint.csv")


save(merged_drb_impute,
     proteinGroups_Zymo_Eclipse45min,
     peptides_Zymo_Eclipse45min,
     file= "P:/EAT_20190926_Zymomona/Zymo_TimeCourse/LC_RawFiles/Oct2020_Technical_Eclispe/txt_DDA_45minHighResMS2_120TimePoint/EATsimpleStats_DRBimputed_50percentfiltered_proteinGroups_peptides_DDA_45minHighResMS2_120TimePoint.RData")
