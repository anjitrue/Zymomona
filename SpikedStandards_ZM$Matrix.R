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

 df <- SpikeMatrixLadder_DIF_AIEIVDQALDRlight
 target_ions <- c("y7", "y8", "y9")

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
  
  for(i in 1:length(target_ions)){
      df_totalIons_sum <- data.frame(df_totalIons_sum, apply(df[which(df$Mass.Feature %in% target_ions[i]), -c(1:2)],2,sum))
  }
  
  df_totalIons_sum_ordered <- df_totalIons_sum[order(df_totalIons_sum$Concentration),]
  colnames(df_totalIons_sum_ordered) <- c(colnames(Sample_meta), "All", "Select","y7", "y8", "y9")
  
  
  
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


#### Aggregate funtions ####

df<- Summed_SpikeMatrixLadder_ETDIGVTGGGQGKheavy
WhichIons <- "y7"

aggregate_ions <- function(df, WhichIons){
  
  colnames(df) <- c(colnames(df)[1:3], "All", "Select","y7", "y8", "y9")
  
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

### Function using aggreagate with dataframe that contains more than 3 replicate injections per concentration

# df<- Summed_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGKheavy
# WhichIons <- "Select"

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

## Merge dataframes with summation of all ions and select ions

# df_all <- All_SpikeMatrixLadder_ETDIGVTGGGQGKlight
# df_select <- Select_SpikeMatrixLadder_ETDIGVTGGGQGKlight

merge_spiked <- function(df_all, df_select, df_y7, df_y8, df_y9){
  df_all$IonsMerged <- rep("All", nrow(df_all))
  df_select$IonsMerged <- rep("Select", nrow(df_all))
  df_y7$IonsMerged <- rep("y7", nrow(df_all))
  df_y8$IonsMerged <- rep("y8", nrow(df_all))
  df_y9$IonsMerged <- rep("y9", nrow(df_all))
  
  df_merged <- rbind(df_all, df_select, df_y7, df_y8, df_y9)
}



#### Plotting light and heavy peptides separately #####

# x <- Summed_SpikeMatrixLadder_ETDIGVTGGGQGKheavy
# All <- "All"
# Select <- "Select"
# y7 <- "y7"
# y8 <- "y8"
# y9 <- "y9"

collated_spikedMatrix <- function(x,All, Select, y7, y8, y9){
  all_x <- aggregate_ions(x, All)
  select_x <- aggregate_ions(x, Select)
  y7_x <- aggregate_ions(x, y7)
  y8_x <- aggregate_ions(x, y8)
  y9_x <- aggregate_ions(x, y9)
  
  merged_spike_Df <- merge_spiked(all_x, select_x, y7_x, y8_x, y9_x)
  
  return(merged_spike_Df)
}

# ETDIGVTGGGQGK native
Merged_SpikeMatrixLadder_ETDIGVTGGGQGKlight <- collated_spikedMatrix(Summed_SpikeMatrixLadder_ETDIGVTGGGQGKlight, 
                                                                     "All", "Select", "y7", "y8", "y9")

type.colors = as.numeric(factor(Merged_SpikeMatrixLadder_ETDIGVTGGGQGKlight$IonsMerged))
ggplot(Merged_SpikeMatrixLadder_ETDIGVTGGGQGKlight, aes(Concentration, Average_Sum, color=IonsMerged)) +
  geom_point()+
  geom_line()+
  theme_light() + 
  labs(title = "Spiked Matrix Curve - ETDIGVTGGGQGK native", subtitle = "Average Intensity Across sample types in triplicate",
       y = "Summed Transitioin Intensity")

# ETDIGVTGGGQGK  standard
Merged_SpikeMatrixLadder_ETDIGVTGGGQGKheavy <-collated_spikedMatrix(Summed_SpikeMatrixLadder_ETDIGVTGGGQGKheavy, 
                                                                    "All", "Select", "y7", "y8", "y9")
type.colors = as.numeric(factor(Merged_SpikeMatrixLadder_ETDIGVTGGGQGKheavy$IonsMerged))
ggplot(Merged_SpikeMatrixLadder_ETDIGVTGGGQGKheavy, aes(Concentration, Average_Sum, color=IonsMerged)) +
  geom_point()+
  geom_line()+
  theme_light() + 
  labs(title = "Spiked Matrix Curve - ETDIGVTGGGQGK standard", subtitle = "Average Intensity Across sample types in triplicate",
       y = "Summed Transitioin Intensity")

type.colors = as.numeric(factor(Merged_SpikeMatrixLadder_ETDIGVTGGGQGKheavy$IonsMerged))
ggplot(Merged_SpikeMatrixLadder_ETDIGVTGGGQGKheavy[-which(Merged_SpikeMatrixLadder_ETDIGVTGGGQGKheavy$IonsMerged == "All"),], aes(Concentration, Average_Sum, color=IonsMerged)) +
  geom_point()+
  geom_line()+
  theme_light() + 
  labs(title = "Spiked Matrix Curve - ETDIGVTGGGQGK standard", subtitle = "Average Inttensity Across sample types in triplicate",
       y = "Summed Transitioin Intensity")

# AIEIVDQALDR  native
Merged_SpikeMatrixLadder_AIEIVDQALDRlight <- collated_spikedMatrix(Summed_SpikeMatrixLadder_AIEIVDQALDRlight, 
                                                                   "All", "Select", "y7", "y8", "y9")

type.colors = as.numeric(factor(Merged_SpikeMatrixLadder_AIEIVDQALDRlight$IonsMerged))
ggplot(Merged_SpikeMatrixLadder_AIEIVDQALDRlight, aes(Concentration, Average_Sum, color=IonsMerged)) +
  geom_point()+
  geom_line()+
  theme_light() + 
  labs(title = "Spiked Matrix Curve - AIEIVDQALDR native", subtitle = "Average Inttensity Across sample types in triplicate",
       y = "Summed Transitioin Intensity")

type.colors = as.numeric(factor(Merged_SpikeMatrixLadder_AIEIVDQALDRlight$IonsMerged))
ggplot(Merged_SpikeMatrixLadder_AIEIVDQALDRlight[-which(Merged_SpikeMatrixLadder_AIEIVDQALDRlight$IonsMerged == "All"),], aes(Concentration, Average_Sum, color=IonsMerged)) +
  geom_point()+
  geom_line()+
  theme_light() + 
  labs(title = "Spiked Matrix Curve - AIEIVDQALDR native", subtitle = "Average Inttensity Across sample types in triplicate",
       y = "Summed Transitioin Intensity")

# AIEIVDQALDR standard
Merged_SpikeMatrixLadder_AIEIVDQALDRheavy <- collated_spikedMatrix(Summed_SpikeMatrixLadder_AIEIVDQALDRheavy, 
                                                                   "All", "Select", "y7", "y8", "y9")

type.colors = as.numeric(factor(Merged_SpikeMatrixLadder_AIEIVDQALDRheavy$IonsMerged))
ggplot(Merged_SpikeMatrixLadder_AIEIVDQALDRheavy, aes(Concentration, Average_Sum, color=IonsMerged)) +
  geom_point()+
  geom_line()+
  theme_light() + 
  labs(title = "Spiked Matrix Curve - AIEIVDQALDR standard", subtitle = "Average Inttensity Across sample types in triplicate",
       y = "Summed Transitioin Intensity")

ggplot(Merged_SpikeMatrixLadder_AIEIVDQALDRheavy[-which(Merged_SpikeMatrixLadder_AIEIVDQALDRheavy$IonsMerged == "All"),], aes(Concentration, Average_Sum, color=IonsMerged)) +
  geom_point()+
  geom_line()+
  theme_light() + 
  labs(title = "Spiked Matrix Curve - AIEIVDQALDR standard", subtitle = "Average Inttensity Across sample types in triplicate",
       y = "Summed Transitioin Intensity")

# FAIMS ETDIGVTGGGQGK  native
Merged_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGKlight <- collated_spikedMatrix(Summed_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGKlight, 
                                                                     "All", "Select", "y7", "y8", "y9")


type.colors = as.numeric(factor(Merged_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGKlight$IonsMerged))
ggplot(Merged_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGKlight, aes(Concentration, Average_Sum, color=IonsMerged)) +
  geom_point()+
  geom_line()+
  theme_light() + 
  labs(title = "Spiked Matrix Curve - ETDIGVTGGGQGK native", subtitle = "Average Intensity Across sample types in triplicate",
       y = "Summed Transitioin Intensity")

# FAIMS ETDIGVTGGGQGK  standard

Merged_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGKheavy <- collated_spikedMatrix(Summed_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGKheavy, 
                                                                           "All", "Select", "y7", "y8", "y9")

type.colors = as.numeric(factor(Merged_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGKheavy$IonsMerged))
ggplot(Merged_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGKheavy, aes(Concentration, Average_Sum, color=IonsMerged)) +
  geom_point()+
  geom_line()+
  theme_light() + 
  labs(title = "Spiked Matrix Curve - FAIMS ETDIGVTGGGQGK standard", subtitle = "Average Intensity Across sample types in triplicate",
       y = "Summed Transitioin Intensity")

ggplot(Merged_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGKheavy[-which(Merged_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGKheavy$IonsMerged == "All"),],
       aes(Concentration, Average_Sum, color=IonsMerged)) +
  geom_point()+
  geom_line()+
  theme_light() + 
  labs(title = "Spiked Matrix Curve - FAIMS ETDIGVTGGGQGK standard", subtitle = "Average Intensity Across sample types in triplicate",
       y = "Summed Transitioin Intensity")

# FAIMS AIEIVDQALDR  native
Merged_FAIMS_SpikeMatrixLadder_AIEIVDQALDRlight <- collated_spikedMatrix(Summed_FAIMS_SpikeMatrixLadder_AIEIVDQALDRlight, 
                                                                         "All", "Select", "y7", "y8", "y9")

type.colors = as.numeric(factor(Merged_FAIMS_SpikeMatrixLadder_AIEIVDQALDRlight$IonsMerged))
ggplot(Merged_FAIMS_SpikeMatrixLadder_AIEIVDQALDRlight, aes(Concentration, Average_Sum, color=IonsMerged)) +
  geom_point()+
  geom_line()+
  theme_light() + 
  labs(title = "Spiked Matrix Curve - FAIMS AIEIVDQALDR native", subtitle = "Average Inttensity Across sample types in triplicate",
       y = "Summed Transitioin Intensity")

type.colors = as.numeric(factor(Merged_FAIMS_SpikeMatrixLadder_AIEIVDQALDRlight$IonsMerged))
ggplot(Merged_FAIMS_SpikeMatrixLadder_AIEIVDQALDRlight[-which(Merged_SpikeMatrixLadder_AIEIVDQALDRlight$IonsMerged == "All"),], aes(Concentration, Average_Sum, color=IonsMerged)) +
  geom_point()+
  geom_line()+
  theme_light() + 
  labs(title = "Spiked Matrix Curve - FAIMS AIEIVDQALDR native", subtitle = "Average Inttensity Across sample types in triplicate",
       y = "Summed Transitioin Intensity")

Merged_FAIMS_SpikeMatrixLadder_AIEIVDQALDRheavy <- collated_spikedMatrix(Summed_FAIMS_SpikeMatrixLadder_AIEIVDQALDRheavy, 
                                                                         "All", "Select", "y7", "y8", "y9")

type.colors = as.numeric(factor(Merged_FAIMS_SpikeMatrixLadder_AIEIVDQALDRheavy$IonsMerged))
ggplot(Merged_FAIMS_SpikeMatrixLadder_AIEIVDQALDRheavy, aes(Concentration, Average_Sum, color=IonsMerged)) +
  geom_point()+
  geom_line()+
  theme_light() + 
  labs(title = "Spiked Matrix Curve - FAIMS AIEIVDQALDR standard", subtitle = "Average Inttensity Across sample types in triplicate",
       y = "Summed Transitioin Intensity")




# df_light <- Summed_SpikeMatrixLadder_ETDIGVTGGGQGKlight
# df_heavy <- Summed_SpikeMatrixLadder_ETDIGVTGGGQGKheavy
# WhichIons <- c("All","Select","y7","y8","y9")

ratio_aggregate <- function(df_light, df_heavy,WhichIons){
  
  # colnames(df_light) <- c(colnames(df_light[,1:3]), WhichIons, "Concentration")
  # colnames(df_heavy) <- c(colnames(df_heavy[,1:3]), WhichIons,"Concentration")
  df_heavy_ions <- which(colnames(df_heavy) %in% WhichIons)
  df_light_ions <- which(colnames(df_light) %in% WhichIons)
  
  ratio <- df_heavy[df_heavy_ions]/df_light[df_light_ions]
  
  df_ratio <- data.frame(df_light$Concentration, ratio)
  
  aggregate_sum_ions <- aggregate(df_ratio[,-1],
                                  list(df_ratio$df_light.Concentration), mean)
  
  colnames(aggregate_sum_ions) <- c("Concentration", "All", "Select", "y7", "y8", "y9")
  aggregate_sum_ions$Aggregate <- rep("Average", nrow(aggregate_sum_ions))
  
  aggregate_std_ions <- aggregate(df_ratio[,-1],
                                        list(df_ratio$df_light.Concentration), sd)
  colnames(aggregate_std_ions) <- c("Concentration", "All", "Select", "y7", "y8", "y9")
  aggregate_std_ions$Aggregate <- rep("Stdev", nrow(aggregate_std_ions))
  
  aggregate_CV_ions <- data.frame(aggregate_std_ions$Concentration ,(aggregate_std_ions[,-c(1,7)]/aggregate_sum_ions[,-c(1,7)])*100)
  colnames(aggregate_CV_ions) <- c( "Concentration", "All", "Select", "y7", "y8", "y9")
  aggregate_CV_ions$Aggregate <- rep("CV", nrow(aggregate_std_ions))
  
  aggregate_ratio_ions <- rbind(aggregate_sum_ions, aggregate_std_ions, aggregate_CV_ions)
  
  return(aggregate_ratio_ions)
}

# df_light <- Summed_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGKlight
# df_heavy <- Summed_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGKheavy


Repliciate_ratio_aggregate <- function(df_light, df_heavy){
  
  colnames(df_light) <- c(colnames(df_light[,1:3]), "All", "Select", "y7", "y8", "y9", "RepConcentration")
  df_light <- df_light[-grep("reinject", df_light$MultiPlexReplicate),]
  colnames(df_heavy) <- c(colnames(df_heavy[,1:3]), "All", "Select", "y7", "y8", "y9", "RepConcentration")
  df_heavy <- df_heavy[-grep("reinject", df_heavy$MultiPlexReplicate),]
  
  ratio <- df_heavy[,-c(1:3,9)]/df_light[,-c(1:3,9)]
  
  df_ratio <- data.frame(df_light$Concentration, ratio)
  
  aggregate_sum_ions <- aggregate(df_ratio[,-1],
                                  list(df_ratio$df_light.Concentration), mean)
  
  colnames(aggregate_sum_ions) <- c("Concentration", "All", "Select", "y7", "y8", "y9")
  aggregate_sum_ions$Aggregate <- rep("Average", nrow(aggregate_sum_ions))
  
  aggregate_std_ions <- aggregate(df_ratio[,-1],
                                  list(df_ratio$df_light.Concentration), sd)
  colnames(aggregate_std_ions) <- c("Concentration", "All", "Select", "y7", "y8", "y9")
  aggregate_std_ions$Aggregate <- rep("Stdev", nrow(aggregate_std_ions))
  
  aggregate_CV_ions <- data.frame(aggregate_std_ions$Concentration ,(aggregate_std_ions[,-c(1,7)]/aggregate_sum_ions[,-c(1,7)])*100)
  colnames(aggregate_CV_ions) <- c( "Concentration", "All", "Select", "y7", "y8", "y9")
  aggregate_CV_ions$Aggregate <- rep("CV", nrow(aggregate_std_ions))
  
  aggregate_ratio_ions <- rbind(aggregate_sum_ions, aggregate_std_ions, aggregate_CV_ions)
  
  return(aggregate_ratio_ions)
}


# No FAIMS
Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK <- ratio_aggregate(Summed_SpikeMatrixLadder_ETDIGVTGGGQGKlight,Summed_SpikeMatrixLadder_ETDIGVTGGGQGKheavy, c("All","Select","y7","y8","y9"))

Ratio_SpikeMatrixLadder_AIEIVDQALDR <- ratio_aggregate(Summed_SpikeMatrixLadder_AIEIVDQALDRlight,Summed_SpikeMatrixLadder_AIEIVDQALDRheavy, c("All","Select","y7","y8","y9"))


# FAIMS
Ratio_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGK <- Repliciate_ratio_aggregate(Summed_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGKlight,Summed_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGKheavy)

#Ratio_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGK_All_ordered <- Ratio_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGK_All %>%
#  mutate(Group.1 = factor(Group.1, level = Ratio_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGK_All$Group.1))



Ratio_FAIMS_SpikeMatrixLadder_AIEIVDQALDR <- Repliciate_ratio_aggregate(Summed_FAIMS_SpikeMatrixLadder_AIEIVDQALDRlight,Summed_FAIMS_SpikeMatrixLadder_AIEIVDQALDRheavy)




#Plot ETDIGVTGGGQGK ratio with and without FAIMMS
Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK$IonMobility <- rep("noFAIMS", nrow(Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK))
Ratio_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGK$IonMobility <- rep("FAIMS", nrow(Ratio_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGK))

#Bind FAIMS and NoFAIMS data
Combined_Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK <- rbind(Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK,Ratio_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGK)
Combined_Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK$Concentration <- c(Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK$Group.1,Ratio_FAIMS_SpikeMatrixLadder_ETDIGVTGGGQGK$Concentration)
Combined_Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK <- Combined_Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK[order(Combined_Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK$Concentration),]

Combined_Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK$Type <- paste0(Combined_Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK$Concentration, Combined_Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK$IonMobility)

Combined_Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK_Avg <- Combined_Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK[grep("Average",
                                                                                                          Combined_Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK$Aggregate),]
Combined_Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK_Stdev <- Combined_Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK[grep("Stdev",
                                                                                                          Combined_Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK$Aggregate),]
Combined_Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK_Avg_long <- melt(Combined_Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK_Avg[,-c(7)], id.vars = c("Concentration", "IonMobility", "Type"))
colnames(Combined_Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK_Avg_long) <- c("Concentration", "IonMobility", "Type",  "variable", "Average_value")
Combined_Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK_Stdev_long <- melt(Combined_Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK_Stdev[,-c(7)], id.vars = c("Concentration", "IonMobility", "Type"))
Combined_Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK_Avg_Stdev_long <- Combined_Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK_Avg_long
Combined_Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK_Avg_Stdev_long$Stdev <- Combined_Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK_Stdev_long$value
#fitlm_ETDIGVTGGGQGK_noFAIMS <- lm(value ~ Concentration, data = Combined_Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK_Avg_long[grep("noFAIMS",Combined_Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK_Avg_long$Type),])
#fitlm_ETDIGVTGGGQGK_FAIMS <- lm(value ~ Concentration, data = Combined_Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK_Avg_long[-grep("noFAIMS",Combined_Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK_Avg_long$Type),])


#predict_ETDIGVTGGGQK <- as.data.frame(cbind(c(predict(fitlm_ETDIGVTGGGQGK_noFAIMS), predict(fitlm_ETDIGVTGGGQGK_FAIMS)), predict_ETDIGVTGGGQK$IonMobility <- c(rep("noFAIMS",6), rep("FAIMS", 6))))
#colnames(predict_ETDIGVTGGGQK) <- c("Predicted", "IonMobility")
#predict_ETDIGVTGGGQK$Predicted <- as.numeric(predict_ETDIGVTGGGQK$Predicted)
#predict_ETDIGVTGGGQK$Concentration <- rep(c(15.625,31.250, 62.5, 125.0, 250.0, 500.0),2)



summary(lm(fitlm_ETDIGVTGGGQGK))$r.squared


Combined_Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK_Avg_Stdev_long$IonMobility <- as.factor(Combined_Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK_Avg_Stdev_long$IonMobility)
PE <- ggplot(Combined_Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK_Avg_Stdev_long[which(Combined_Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK_Avg_Stdev_long$IonMobility == "FAIMS"),],
             aes(x= Concentration, y=Average_value, group = variable)) +
  geom_point(aes(color = variable, shape = variable), size = 6, alpha = 0.7) +
  scale_shape_manual(values = 1:nlevels(Combined_Ratio_SpikeMatrixLadder_ETDIGVTGGGQGK_Avg_Stdev_long$variable))+
  #geom_line(aes(color = variable)) +
  geom_errorbar(aes(ymin = Average_value-Stdev, ymax = Average_value+Stdev),
                alpha = 0.2,
                width = 0,
                size= 10,
                color= "#6D696F") +
  #geom_line(data=predict_ETDIGVTGGGQK, aes(x=Concentration,y=Predicted, color = IonMobility))+
  theme_classic()+
  scale_y_continuous(breaks = seq(0,0.6,0.1), limit = c(0,0.6))+
  #theme(axis.text.x = element_text(angle = 0, hjust = 1))+
  labs(title = "Standard Spiked in Matrix", subtitle = "Heavy to Light Ratio - ETDIGVTGGGQGK",
       x = "Concentration (fmol/uL)", y = "heavy to light Ratio")



#Plot AIEIVDQALDR ratio with and without FAIMS
Ratio_SpikeMatrixLadder_AIEIVDQALDR$IonMobility <- rep("noFAIMS", nrow(Ratio_SpikeMatrixLadder_AIEIVDQALDR))
Ratio_FAIMS_SpikeMatrixLadder_AIEIVDQALDR$IonMobility <- rep("FAIMS", nrow(Ratio_FAIMS_SpikeMatrixLadder_AIEIVDQALDR))

Combined_Ratio_SpikeMatrixLadder_AIEIVDQALDR <- rbind(Ratio_SpikeMatrixLadder_AIEIVDQALDR,Ratio_FAIMS_SpikeMatrixLadder_AIEIVDQALDR)
Combined_Ratio_SpikeMatrixLadder_AIEIVDQALDR$Concentration <- c(Ratio_SpikeMatrixLadder_AIEIVDQALDR$Group.1,Ratio_FAIMS_SpikeMatrixLadder_AIEIVDQALDR$Concentration)
Combined_Ratio_SpikeMatrixLadder_AIEIVDQALDR <- Combined_Ratio_SpikeMatrixLadder_AIEIVDQALDR[order(Combined_Ratio_SpikeMatrixLadder_AIEIVDQALDR$Concentration),]

Combined_Ratio_SpikeMatrixLadder_AIEIVDQALDR$Type <- paste0(Combined_Ratio_SpikeMatrixLadder_AIEIVDQALDR$Concentration, Combined_Ratio_SpikeMatrixLadder_AIEIVDQALDR$IonMobility)

Combined_Ratio_SpikeMatrixLadder_AIEIVDQALDR_Avg <- Combined_Ratio_SpikeMatrixLadder_AIEIVDQALDR[grep("Average",
                                                                                                          Combined_Ratio_SpikeMatrixLadder_AIEIVDQALDR$Aggregate),]

Combined_Ratio_SpikeMatrixLadder_AIEIVDQALDR_Stdev <- Combined_Ratio_SpikeMatrixLadder_AIEIVDQALDR[grep("Stdev",
                                                                                                            Combined_Ratio_SpikeMatrixLadder_AIEIVDQALDR$Aggregate),]
Combined_Ratio_SpikeMatrixLadder_AIEIVDQALDR_Avg_long <- melt(Combined_Ratio_SpikeMatrixLadder_AIEIVDQALDR_Avg[,-c(7)], id.vars = c("Concentration", "IonMobility", "Type"))
colnames(Combined_Ratio_SpikeMatrixLadder_AIEIVDQALDR_Avg_long) <- c("Concentration", "IonMobility", "Type",  "variable", "Average_value")
Combined_Ratio_SpikeMatrixLadder_AIEIVDQALDR_Stdev_long <- melt(Combined_Ratio_SpikeMatrixLadder_AIEIVDQALDR_Stdev[,-c(7)], id.vars = c("Concentration", "IonMobility", "Type"))
Combined_Ratio_SpikeMatrixLadder_AIEIVDQALDR_Avg_Stdev_long <- Combined_Ratio_SpikeMatrixLadder_AIEIVDQALDR_Avg_long
Combined_Ratio_SpikeMatrixLadder_AIEIVDQALDR_Avg_Stdev_long$Stdev <- Combined_Ratio_SpikeMatrixLadder_AIEIVDQALDR_Stdev_long$value

#fitlm_AIEIVDQALDR_noFAIMS <- lm(value ~ Concentration, data = Combined_Ratio_SpikeMatrixLadder_AIEIVDQALDR_Avg_long[grep("noFAIMS",Combined_Ratio_SpikeMatrixLadder_AIEIVDQALDR_Avg_long$Type),])
#fitlm_AIEIVDQALDR_FAIMS <- lm(value ~ Concentration, data = Combined_Ratio_SpikeMatrixLadder_AIEIVDQALDR_Avg_long[-grep("noFAIMS",Combined_Ratio_SpikeMatrixLadder_AIEIVDQALDR_Avg_long$Type),])


# predict_ETDIGVTGGGQK <- as.data.frame(cbind(c(predict(fitlm_AIEIVDQALDR_noFAIMS), predict(fitlm_AIEIVDQALDR_FAIMS)), predict_AIEIVDQALDR$IonMobility <- c(rep("noFAIMS",6), rep("FAIMS", 6))))
# colnames(predict_AIEIVDQALDR) <- c("Predicted", "IonMobility")
# predict_AIEIVDQALDR$Predicted <- as.numeric(predict_AIEIVDQALDR$Predicted)
# predict_AIEIVDQALDR$Concentration <- rep(c(15.625,31.250, 62.5, 125.0, 250.0, 500.0),2)



#summary(lm(fitlm_ETDIGVTGGGQGK))$r.squared


Combined_Ratio_SpikeMatrixLadder_AIEIVDQALDR_Avg_Stdev_long$IonMobility <- as.factor(Combined_Ratio_SpikeMatrixLadder_AIEIVDQALDR_Avg_Stdev_long$IonMobility)
PA <- ggplot(Combined_Ratio_SpikeMatrixLadder_AIEIVDQALDR_Avg_Stdev_long[which(Combined_Ratio_SpikeMatrixLadder_AIEIVDQALDR_Avg_Stdev_long$IonMobility == "FAIMS"),], 
             aes(x= Concentration, y=Average_value, group = variable)) +
  geom_point(aes(color = variable, shape = variable), size = 6, alpha = 0.7) +
  scale_shape_manual(values = 1:nlevels(Combined_Ratio_SpikeMatrixLadder_AIEIVDQALDR_Avg_Stdev_long$variable))+
  geom_errorbar(aes(ymin = Average_value-Stdev, ymax = Average_value+Stdev),
                alpha = 0.2,
                width = 0,
                size= 10,
                color= "#6D696F") +
  #geom_line(data=predict_ETDIGVTGGGQK, aes(x=Concentration,y=Predicted, color = IonMobility))+
  scale_y_continuous(breaks = seq(0,18,2), limit = c(0,18))+
  theme_classic()+
  #theme(axis.text.x = element_text(angle = 0, hjust = 1))+
  labs(title = "Standard Spiked in Matrix", subtitle = "Heavy to Light Ratio - AIEIVDQALDR",
       x = "Concentration (fmol/uL)", y = "heavy to light Ratio")

PEgrob <- ggplotGrob(PE)
PAgrob <- ggplotGrob(PA)

pdf("F:/Projects/Proteomics/Zymomona/FAIMS/Figures/FromR/SpectraComparisons_PROSIT_HIGHRES_LOWRES_IW6_withErrorBars.pdf", width = 20, height = 10)
grid.arrange(PEgrob, PAgrob, ncol =2)
dev.off()
