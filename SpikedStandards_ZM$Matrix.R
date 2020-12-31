#### Load Data ####

SpikeMatrixLadder_DIF_ETDIGVTGGGQGKlight <- read.csv("P:/EAT_20190926_Zymomona/Zymo_TimeCourse/Infusions/SpikedMatrixLadder_DirectInfusion/MultiPlexReplicates/CompiledFragmentResults_ETDIGVTGGGQGK_SpikedMatrix_Replicates_IsotopeDifferenceRecalculated.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

SpikeMatrixLadder_DIF_ETDIGVTGGGQGKheavy <- read.csv("P:/EAT_20190926_Zymomona/Zymo_TimeCourse/Infusions/SpikedMatrixLadder_DirectInfusion/MultiPlexReplicates/CompiledFragmentResults_ETDIGVTGGGQGKheavy_SpikedMatrix_Replicates_IsotopeDifferenceRecalculated.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

SpikeMatrixLadder_DIF_AIEIVDQALDRlight <- read.csv("P:/EAT_20190926_Zymomona/Zymo_TimeCourse/Infusions/SpikedMatrixLadder_DirectInfusion/MultiPlexReplicates/CompiledFragmentResults_AIEIVDQALDR_DirectInfusion_SpikedMatrix_Replicates_IsotopeDifferenceRecalculated.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

SpikeMatrixLadder_DIF_AIEIVDQALDRheavy <- read.csv("P:/EAT_20190926_Zymomona/Zymo_TimeCourse/Infusions/SpikedMatrixLadder_DirectInfusion/MultiPlexReplicates/CompiledFragmentResults_AIEIVDQALDRheavy_DirectInfusion_SpikedMatrix_Replicates_IsotopeDifferenceRecalculated.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

### FAIMS 

SpikeMatrixLadder_FAIMS_ETDIGVTGGGQGKlight <- read.csv("P:/EAT_20190926_Zymomona/Zymo_TimeCourse/Infusions/SpikedMatrixLadder_FAIMS/MultiPlexReplicates/CompiledFragmentResults_ETDIGVTGGGQGK_FAIMS_SpikedMatrix_Replicates_IsotopeDifferenceRecalculated.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

SpikeMatrixLadder_FAIMS_ETDIGVTGGGQGKheavy <- read.csv("P:/EAT_20190926_Zymomona/Zymo_TimeCourse/Infusions/SpikedMatrixLadder_FAIMS/MultiPlexReplicates/CompiledFragmentResults_ETDIGVTGGGQGKheavy_FAIMS_SpikedMatrix_Replicates_IsotopeDifferenceRecalculated.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

SpikeMatrixLadder_FAIMS_AIEIVDQALDRlight <- read.csv("P:/EAT_20190926_Zymomona/Zymo_TimeCourse/Infusions/SpikedMatrixLadder_FAIMS/MultiPlexReplicates/CompiledFragmentResults_AIEIVDQALDR_FAIMS_DirectInfusion_SpikedMatrix_Replicates_IsotopeDifferenceRecalculated.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

SpikeMatrixLadder_FAIMS_AIEIVDQALDRheavy <- read.csv("P:/EAT_20190926_Zymomona/Zymo_TimeCourse/Infusions/SpikedMatrixLadder_FAIMS/MultiPlexReplicates/CompiledFragmentResults_AIEIVDQALDRheavy_FAIMS_DirectInfusion_SpikedMatrix_Replicates_IsotopeDifferenceRecalculated.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)


#### SpikedMatriix_Formatting #####

# df <- SpikeMatrixLadder_FAIMS_ETDIGVTGGGQGKheavy
# target_ions <- c("y7", "y8", "y9")

Induction_Summed_precusor <- function(df, target_ions){
  
  #Extract sample names and reformat
  RAW_Sample_names <- colnames(df[,-c(1:2)])
  
  for(i in 1:length(RAW_Sample_names)){
    
    if(length(grep("2020", RAW_Sample_names[i])) ==!0){
      replicate_number <- sub(".*Multiplex","",RAW_Sample_names[i])
      replicate_number <- sub("_.*", "", replicate_number)
      reinject_name <- paste0("Multiplex","reinject", replicate_number)
      RAW_Sample_names[i] <- sub("Multiplex.*",reinject_name, RAW_Sample_names[i])
      
    }
  }
  
  Concentration <- sub("_spiked.*", "", RAW_Sample_names)
  Concentration <- sub("X", "", Concentration)
  Concentration <- sub("_FAIMS", "", Concentration)
  #Concentration <- gsub('.{1}$', '', Concentration)
  Concentration <- as.numeric(sub("_", ".", Concentration))
  
  MultiPlexReplicate <- sub(".*Matrix_", "", RAW_Sample_names)
  MultiPlexReplicate <- sub(".raw", "", MultiPlexReplicate)
  MultiPlexReplicate <- sub("Multiplex","", MultiPlexReplicate)
  
  Sample_meta <- data.frame(RAW_Sample_names, Concentration, MultiPlexReplicate)
  
  df_totalIons_sum <- data.frame(Sample_meta,apply(df[,-c(1:2)],2,sum),
                                 apply(df[which(df$Mass.Feature %in% target_ions),-c(1:2)],2,sum))
  
  df_totalIons_sum_ordered <- df_totalIons_sum[order(df_totalIons_sum$Concentration),]
  colnames(df_totalIons_sum_ordered) <- c(colnames(Sample_meta), "Sum_All_Ions", "Sum_SelectIons")
  
  # df_totalIons_sum_ordered <- df_totalIons_sum %>%
  # mutate(Concentttration = factor(Concentration, level = unique(df_totalIons_sum$Concentration)))
  return(df_totalIons_sum_ordered)
  
}

Summed_SpikeMatrixLadder_ETDIGVTGGGQGKlight <-  Induction_Summed_precusor(SpikeMatrixLadder_DIF_ETDIGVTGGGQGKlight,c("y7", "y8", "y9"))

Summed_SpikeMatrixLadder_ETDIGVTGGGQGKheavy <-  Induction_Summed_precusor(SpikeMatrixLadder_DIF_ETDIGVTGGGQGKheavy,c("y7", "y8", "y9"))

Summed_SpikeMatrixLadder_AIEIVDQALDRlight <-  Induction_Summed_precusor(SpikeMatrixLadder_DIF_AIEIVDQALDRlight,c("y7", "y8", "y9"))

Summed_SpikeMatrixLadder_AIEIVDQALDRheavy <-  Induction_Summed_precusor(SpikeMatrixLadder_DIF_AIEIVDQALDRheavy,c("y7", "y8", "y9"))

### FAIMS
Summed_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGKlight <-  Induction_Summed_precusor(SpikeMatrixLadder_FAIMS_ETDIGVTGGGQGKlight ,c("y7", "y8", "y9"))
Summed_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGKlight$RepConcentration <- paste0(Summed_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGKlight$Concentration,gsub('.{1}$','',Summed_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGKlight$MultiPlexReplicate))

Summed_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGKheavy <-  Induction_Summed_precusor(SpikeMatrixLadder_FAIMS_ETDIGVTGGGQGKheavy,c("y7", "y8", "y9"))
Summed_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGKheavy$RepConcentration <- paste0(Summed_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGKheavy$Concentration,gsub('.{1}$','',Summed_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGKheavy$MultiPlexReplicate))

Summed_FAIMS_SpikeMatrixLadder_AIEIVDQALDRlight <-  Induction_Summed_precusor(SpikeMatrixLadder_FAIMS_AIEIVDQALDRlight ,c("y7", "y8", "y9"))
Summed_FAIMS_SpikeMatrixLadder_AIEIVDQALDRlight$RepConcentration <- paste0(Summed_FAIMS_SpikeMatrixLadder_AIEIVDQALDRlight$Concentration,gsub('.{1}$','',Summed_FAIMS_SpikeMatrixLadder_AIEIVDQALDRlight$MultiPlexReplicate))

Summed_FAIMS_SpikeMatrixLadder_AIEIVDQALDRheavy <-  Induction_Summed_precusor(SpikeMatrixLadder_FAIMS_AIEIVDQALDRheavy,c("y7", "y8", "y9"))
Summed_FAIMS_SpikeMatrixLadder_AIEIVDQALDRheavy$RepConcentration <- paste0(Summed_FAIMS_SpikeMatrixLadder_AIEIVDQALDRheavy$Concentration,gsub('.{1}$','',Summed_FAIMS_SpikeMatrixLadder_AIEIVDQALDRheavy$MultiPlexReplicate))



```{r,echo=FALSE}

df<- Summed_SpikeMatrixLadder_ETDIGVTGGGQGKheavy
WhichIons <- "Select"

aggregate_ions <- function(df, WhichIons){
  
  colnames(df) <- c(colnames(df)[1:3], "All", "Select")
  
  aggregate_sum_ions <- aggregate(df[,which(colnames(df) == WhichIons)],
                                  list(df$Concentration), mean)
  
  aggregate_sum_ions$stdev <- aggregate(df[,which(colnames(df) == WhichIons)],
                                        list(df$Concentration), sd)[[2]]
  
  aggregate_sum_ions$log2Average <- log2(aggregate_sum_ions$x)
  
  aggregate_sum_ions$log2Stdev <- aggregate(df[,which(colnames(df) == WhichIons)],
                                            list(df$Concentration), function(x) sd(log2(x)))[[2]]
  aggregate_sum_ions$CV <- aggregate_sum_ions$stdev/aggregate_sum_ions$x*100
  
  colnames(aggregate_sum_ions) <- c("Concentration", "Average_Sum", "Stdev_Sum", "Log2Average_Sum", "Stdev_Log2Sum", "%CV")
  
  # ggplot(aggregate_sum_ions, aes(Concentration, Log2Average_Sum))+
  # geom_point()
  
  return(aggregate_sum_ions)
}

df<- Summed_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGKheavy
WhichIons <- "Select"

Replicate_injection_aggregate_ions <- function(df, WhichIons){
  
  colnames(df) <- c(colnames(df)[1:3], "All", "Select","RepConcentration")
  
  aggregate_sum_ions <- aggregate(df[,which(colnames(df) == WhichIons)],
                                  list(df$RepConcentration), mean)
  
  aggregate_sum_ions$stdev <- aggregate(df[,which(colnames(df) == WhichIons)],
                                        list(df$RepConcentration), sd)[[2]]
  
  aggregate_sum_ions$log2Average <- log2(aggregate_sum_ions$x)
  
  aggregate_sum_ions$log2Stdev <- aggregate(df[,which(colnames(df) == WhichIons)],
                                            list(df$RepConcentration), function(x) sd(log2(x)))[[2]]
  aggregate_sum_ions$CV <- aggregate_sum_ions$stdev/aggregate_sum_ions$x*100
  
  colnames(aggregate_sum_ions) <- c("RepConcentration", "Average_Sum", "Stdev_Sum", "Log2Average_Sum", "Stdev_Log2Sum", "%CV")
  
  aggregate_sum_ions$Concentration <- as.numeric(sub("reinject", "", aggregate_sum_ions$RepConcentration))
  
  aggregate_sum_ions_ordered <- aggregate_sum_ions[order(aggregate_sum_ions$Concentration),]
  
  # ggplot(aggregate_sum_ions, aes(Concentration, Log2Average_Sum))+
  # geom_point()
  
  return(aggregate_sum_ions_ordered)
}

df_all <- All_SpikeMatrixLadder_ETDIGVTGGGQGKlight
df_select <- Select_SpikeMatrixLadder_ETDIGVTGGGQGKlight

merge_spiked <- function(df_all, df_select){
  df_all$IonsMerged <- rep("All", nrow(df_all))
  df_select$IonsMerged <- rep("Select", nrow(df_all))
  
  df_merged <- rbind(df_all, df_select)
}


```


#### Coefficient of Variation Calculations #####

All_SpikeMatrixLadder_ETDIGVTGGGQGKlight <- aggregate_ions(Summed_SpikeMatrixLadder_ETDIGVTGGGQGKlight, "All")
Select_SpikeMatrixLadder_ETDIGVTGGGQGKlight <- aggregate_ions(Summed_SpikeMatrixLadder_ETDIGVTGGGQGKlight, "Select")
Merged_SpikeMatrixLadder_ETDIGVTGGGQGKlight <- merge_spiked(All_SpikeMatrixLadder_ETDIGVTGGGQGKlight, Select_SpikeMatrixLadder_ETDIGVTGGGQGKlight)


type.colors = as.numeric(factor(Merged_SpikeMatrixLadder_ETDIGVTGGGQGKlight$IonsMerged))
ggplot(Merged_SpikeMatrixLadder_ETDIGVTGGGQGKlight, aes(Concentration, Average_Sum, color=IonsMerged)) +
  geom_point()+
  geom_line()+
  theme_light() + 
  labs(title = "Spiked Matrix Curve - ETDIGVTGGGQGK light", subtitle = "Average Inttensity Across sample types in triplicate",
       y = "Summed Transitioin Intensity")


All_SpikeMatrixLadder_ETDIGVTGGGQGKheavy <- aggregate_ions(Summed_SpikeMatrixLadder_ETDIGVTGGGQGKheavy, "All")
Select_SpikeMatrixLadder_ETDIGVTGGGQGKheavy <- aggregate_ions(Summed_SpikeMatrixLadder_ETDIGVTGGGQGKheavy, "Select")
Merged_SpikeMatrixLadder_ETDIGVTGGGQGKheavy <- merge_spiked(All_SpikeMatrixLadder_ETDIGVTGGGQGKheavy, Select_SpikeMatrixLadder_ETDIGVTGGGQGKheavy)

type.colors = as.numeric(factor(Merged_SpikeMatrixLadder_ETDIGVTGGGQGKheavy$IonsMerged))
ggplot(Merged_SpikeMatrixLadder_ETDIGVTGGGQGKheavy, aes(Concentration, Average_Sum, color=IonsMerged)) +
  geom_point()+
  geom_line()+
  theme_light() + 
  labs(title = "Spiked Matrix Curve - ETDIGVTGGGQGK heavy", subtitle = "Average Inttensity Across sample types in triplicate",
       y = "Summed Transitioin Intensity")

All_SpikeMatrixLadder_AIEIVDQALDRlight <- aggregate_ions(Summed_SpikeMatrixLadder_AIEIVDQALDRlight, "All")
Select_SpikeMatrixLadder_AIEIVDQALDRlight <- aggregate_ions(Summed_SpikeMatrixLadder_AIEIVDQALDRlight, "Select")

All_SpikeMatrixLadder_AIEIVDQALDRheavy <- aggregate_ions(Summed_SpikeMatrixLadder_AIEIVDQALDRheavy, "All")
Select_SpikeMatrixLadder_AIEIVDQALDRheavy <- aggregate_ions(Summed_SpikeMatrixLadder_AIEIVDQALDRheavy, "Select")

#FAIMS

All_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGKlight <- Replicate_injection_aggregate_ions(Summed_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGKlight, "All")
Select_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGKlight <- Replicate_injection_aggregate_ions(Summed_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGKlight, "Select")
Merged_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGKlight <- merge_spiked(All_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGKlight, Select_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGKlight)


type.colors = as.numeric(factor(Merged_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGKlight$IonsMerged))
ggplot(Merged_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGKlight, aes(Concentration, Average_Sum, color=IonsMerged)) +
  geom_point()+
  geom_line()+
  theme_light() + 
  labs(title = "Spiked Matrix Curve - ETDIGVTGGGQGK heavy", subtitle = "Average Inttensity Across sample types in triplicate",
       y = "Summed Transitioin Intensity")

All_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGKheavy <- Replicate_injection_aggregate_ions(Summed_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGKheavy, "All")
Select_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGKheavy <- Replicate_injection_aggregate_ions(Summed_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGKheavy, "Select")
Merged_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGKheavy <- merge_spiked(All_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGKheavy, Select_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGKheavy)


type.colors = as.numeric(factor(Merged_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGKheavy$IonsMerged))
ggplot(Merged_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGKheavy, aes(Concentration, Average_Sum, color=IonsMerged)) +
  geom_point()+
  geom_line()+
  theme_light() + 
  labs(title = "Spiked Matrix Curve - ETDIGVTGGGQGK heavy", subtitle = "Average Inttensity Across sample types in triplicate",
       y = "Summed Transitioin Intensity")

ggplot(Select_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGKheavy, aes(Concentration, Average_Sum)) +
  geom_point()+
  theme_light() + 
  labs(title = "Spiked Matrix Curve - ETDIGVTGGGQGK heavy", subtitle = "Average Inttensity Across sample types in triplicate",
       y = "Summed Transitioin Intensity")

All_FAIMS_SpikeMatrixLadder_AIEIVDQALDRlight <- Replicate_injection_aggregate_ions(Summed_FAIMS_SpikeMatrixLadder_AIEIVDQALDRlight, "All")
Select_SpikeMatrixLadder_AIEIVDQALDRlight <- Replicate_injection_aggregate_ions(Summed_FAIMS_SpikeMatrixLadder_AIEIVDQALDRlight, "Select")

All_FAIMS_SpikeMatrixLadder_AIEIVDQALDRheavy <- Replicate_injection_aggregate_ions(Summed_FAIMS_SpikeMatrixLadder_AIEIVDQALDRheavy, "All")
Select_FAIMS_SpikeMatrixLadder_AIEIVDQALDRheavy <- Replicate_injection_aggregate_ions(Summed_FAIMS_SpikeMatrixLadder_AIEIVDQALDRheavy, "Select")

```

```{r,echo=FALSE}

# df_light <- Summed_SpikeMatrixLadder_ETDIGVTGGGQGKlight
# df_heavy <- Summed_SpikeMatrixLadder_ETDIGVTGGGQGKheavy
# WhichIons <- "All"

ratio_aggregate <- function(df_light, df_heavy, WhichIons){
  
  colnames(df_light) <- c(colnames(df_light[,1:3]), "All", "Select")
  colnames(df_heavy) <- c(colnames(df_heavy[,1:3]), "All", "Select")
  
  ratio <- df_heavy[,which(colnames(df_heavy) == WhichIons)]/df_light[,which(colnames(df_light) == WhichIons)]
  
  df_ratio <- data.frame(df_light$Concentration, ratio)
  
  aggregate_sum_ions <- aggregate(df_ratio$ratio,
                                  list(df_ratio$df_light.Concentration), mean)
  
  aggregate_sum_ions$stdev <- aggregate(df_ratio$ratio,
                                        list(df_ratio$df_light.Concentration), sd)[[2]]
  
  aggregate_sum_ions$CV <- aggregate_sum_ions$stdev/aggregate_sum_ions$x*100
  
  return(aggregate_sum_ions)
}

# df_light <- Summed_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGKlight
# df_heavy <- Summed_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGKheavy
# WhichIons <- "All"

Repliciate_ratio_aggregate <- function(df_light, df_heavy, WhichIons){
  
  colnames(df_light) <- c(colnames(df_light[,1:3]), "All", "Select", "RepConcentration")
  colnames(df_heavy) <- c(colnames(df_heavy[,1:3]), "All", "Select", "RepConcentration")
  
  ratio <- df_heavy[,which(colnames(df_heavy) == WhichIons)]/df_light[,which(colnames(df_light) == WhichIons)]
  
  df_ratio <- data.frame(df_light$Concentration, ratio)
  
  aggregate_sum_ions <- aggregate(df_ratio$ratio,
                                  list(df_ratio$df_light.Concentration), mean)
  
  aggregate_sum_ions$stdev <- aggregate(df_ratio$ratio,
                                        list(df_ratio$df_light.Concentration), sd)[[2]]
  
  aggregate_sum_ions$CV <- aggregate_sum_ions$stdev/aggregate_sum_ions$x*100
  aggregate_sum_ions$Concentration <- as.numeric(sub("reinject", "", aggregate_sum_ions$Group.1))
  
  aggregate_sum_ions_ordered <- aggregate_sum_ions[order(aggregate_sum_ions$Concentration),]
  
  return(aggregate_sum_ions_ordered)
}

```

```{r, echo=FALSE}

Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK_All <- ratio_aggregate(Summed_SpikeMatrixLadder_ETDIGVTGGGQGKlight,Summed_SpikeMatrixLadder_ETDIGVTGGGQGKheavy,"All")

Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK_Select <- ratio_aggregate(Summed_SpikeMatrixLadder_ETDIGVTGGGQGKlight,Summed_SpikeMatrixLadder_ETDIGVTGGGQGKheavy,"Select")


Ratio_SpikeMatrixLadder_AIEIVDQALDR_All <- ratio_aggregate(Summed_SpikeMatrixLadder_AIEIVDQALDRlight,Summed_SpikeMatrixLadder_AIEIVDQALDRheavy,"All")

Ratio_SpikeMatrixLadder_AIEIVDQALDR_Select <- ratio_aggregate(Summed_SpikeMatrixLadder_AIEIVDQALDRlight,Summed_SpikeMatrixLadder_AIEIVDQALDRheavy,"Select")

# FAIMS
Ratio_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGK_All <- Repliciate_ratio_aggregate(Summed_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGKlight,Summed_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGKheavy,"All")

Ratio_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGK_All_ordered <- Ratio_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGK_All %>%
  mutate(Group.1 = factor(Group.1, level = Ratio_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGK_All$Group.1))

Ratio_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGK_Select <- Repliciate_ratio_aggregate(Summed_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGKlight,Summed_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGKheavy,"Select")



Ratio_FAIMS_SpikeMatrixLadder_AIEIVDQALDR_All <- Repliciate_ratio_aggregate(Summed_FAIMS_SpikeMatrixLadder_AIEIVDQALDRlight,Summed_FAIMS_SpikeMatrixLadder_AIEIVDQALDRheavy,"All")

Ratio_FAIMS_SpikeMatrixLadder_AIEIVDQALDR_Select <- Repliciate_ratio_aggregate(Summed_FAIMS_SpikeMatrixLadder_AIEIVDQALDRlight,Summed_FAIMS_SpikeMatrixLadder_AIEIVDQALDRheavy,"Select")

```


```{r, echo=FALSE}
Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK_All$IonMobility <- rep("noFAIMS", nrow(Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK_All))
Ratio_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGK_All$IonMobility <- rep("FAIMS", nrow(Ratio_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGK_All))

Combined_Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK_All <- rbind(Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK_All,Ratio_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGK_All[,-5])
Combined_Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK_All$Concentration <- c(Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK_All$Group.1,Ratio_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGK_All$Concentration)
Combined_Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK_All <- Combined_Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK_All[order(Combined_Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK_All$Concentration),]

Combined_Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK_All$Type <- paste0(Combined_Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK_All$Group.1, Combined_Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK_All$IonMobility)

order_concentration <- factor(Combined_Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK_All$Group.1, levels = c("15.625", "31.25", "62.5", "125", "125reinject", 
                                                                                                     "250", "250reinject", "500", "500reinject"))


Combined_Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK_All_ordered <- Combined_Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK_All %>%
  mutate(Group.1 = factor(Group.1, level = levels(order_concentration)))

fitlm_ETDIGVTGGGQGK_noFAIMS <- lm(x ~ Concentration, data = Combined_Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK_All_ordered[grep("noFAIMS",Combined_Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK_All_ordered$Type),])
fitlm_ETDIGVTGGGQGK_FAIMS <- lm(x ~ Concentration, data = Combined_Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK_All_ordered[-grep("noFAIMS",Combined_Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK_All_ordered$Type),])


predict_ETDIGVTGGGQK <- as.data.frame(cbind(c(predict(fitlm_ETDIGVTGGGQGK_noFAIMS), predict(fitlm_ETDIGVTGGGQGK_FAIMS)), predict_ETDIGVTGGGQK$IonMobility <- c(rep("noFAIMS",6), rep("FAIMS", 6))))
colnames(predict_ETDIGVTGGGQK) <- c("Predicted", "IonMobility")
predict_ETDIGVTGGGQK$Predicted <- as.numeric(predict_ETDIGVTGGGQK$Predicted)
predict_ETDIGVTGGGQK$Concentration <- rep(c(15.625,31.250, 62.5, 125.0, 250.0, 500.0),2)



summary(lm(fitlm_ETDIGVTGGGQGK))$r.squared


Combined_Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK_All_ordered$IonMobility <- as.factor(Combined_Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK_All_ordered$IonMobility)
ggplot(Combined_Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK_All_ordered, aes(x= Concentration, y=x, group = IonMobility)) +
  geom_point(aes(color = IonMobility), size = 4, alpha = 0.5) +
  #geom_line(aes(color = "black")) +
  geom_errorbar(aes(ymin = x-stdev, ymax = x+stdev),
                alpha = 0.2,
                width = 0,
                size= 10,
                color= "#6D696F") +
  geom_line(data=predict_ETDIGVTGGGQK, aes(x=Concentration,y=Predicted, color = IonMobility))+
  theme_classic()+
  #theme(axis.text.x = element_text(angle = 0, hjust = 1))+
  labs(title = "Standard Spiked in Matrix", subtitle = "Heavy to Light Ratio - ETDIGVTGGGQGK",
       x = "Concentration (fmol/uL)", y = "heavy to light Ratio")

```