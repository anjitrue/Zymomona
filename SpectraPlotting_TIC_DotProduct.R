library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(corrplot)
library(ggcorrplot)
library(dplyr)


##### Library Comparison ####
ProteinProspector_AIEIVQALDR <- read.csv("F:/Projects/Proteomics/Zymomona/FAIMS/DataAnalysis/ProteinProspector_AIEIVDQALDR.csv", 
                                         header = TRUE, sep = ",", stringsAsFactors = FALSE)
ProteinProspector_ETDIGVTGGGQGK <- read.csv("F:/Projects/Proteomics/Zymomona/FAIMS/DataAnalysis/ProteinProspector_ETDIGVTGGGQGK.csv", 
                                            header = TRUE, sep = ",", stringsAsFactors = FALSE)

#Library Spectra
#PROSIT
PROSIT_IspH_AIIEIVDALDR_Spectra <- read.csv("F:/Projects/Proteomics/Zymomona/FAIMS/DataAnalysis/LibraryComparison/PROSIT_LibrarySpectra_AIEIVDQALDR.csv", 
                                            header = TRUE, sep = ",", stringsAsFactors = FALSE)
colnames(PROSIT_IspH_AIIEIVDALDR_Spectra) <- c("m.z", "Relative.Abundance")

PROSIT_IspG_ETDIGVTGGGGQGK_Spectra <- read.csv("F:/Projects/Proteomics/Zymomona/FAIMS/DataAnalysis/LibraryComparison/PROSIT_LibrarySpectra_ETDIGVTGGGGQGK.csv", 
                                               header = TRUE, sep = ",", stringsAsFactors = FALSE)
colnames(PROSIT_IspG_ETDIGVTGGGGQGK_Spectra) <- c("m.z", "Relative.Abundance")

#### HIGH_Res data upload ####
High_Res_IspH_AIIEIVDALDR_Spectra <- read.csv("F:/Projects/Proteomics/Zymomona/FAIMS/DataAnalysis/LibraryComparison/HighRes_Intensity_LibrarySpectra_AIEIVDQALDR.csv", 
                                              header = TRUE, sep = ",", stringsAsFactors = FALSE)
High_Res_IspH_AIIEIVDALDR_Spectra.RelativeAbundance <- High_Res_IspH_AIIEIVDALDR_Spectra
High_Res_IspH_AIIEIVDALDR_Spectra.RelativeAbundance$Relative.Abundance <- High_Res_IspH_AIIEIVDALDR_Spectra.RelativeAbundance$Intensity/(max(High_Res_IspH_AIIEIVDALDR_Spectra.RelativeAbundance$Intensity))
colnames(High_Res_IspH_AIIEIVDALDR_Spectra.RelativeAbundance) <- c("m.z", "Intensity", "Relative.Abundance")



High_Res_IspG_ETDIGVTGGGGQGK_Spectra <- read.csv("F:/Projects/Proteomics/Zymomona/FAIMS/DataAnalysis/LibraryComparison/HighRes_LibrarySpectra_ETDIGVTGGGGQGK.csv", 
                                                 header = TRUE, sep = ",", stringsAsFactors = FALSE)
High_Res_IspG_ETDIGVTGGGGQGK_Spectra.RelativeAbundance = High_Res_IspG_ETDIGVTGGGGQGK_Spectra
High_Res_IspG_ETDIGVTGGGGQGK_Spectra.RelativeAbundance$Relative.Abundance <- High_Res_IspG_ETDIGVTGGGGQGK_Spectra.RelativeAbundance$Intensity/(max(High_Res_IspG_ETDIGVTGGGGQGK_Spectra.RelativeAbundance$Intensity))
colnames(High_Res_IspG_ETDIGVTGGGGQGK_Spectra.RelativeAbundance) <- c("m.z", "Intensity", "Relative.Abundance")



#### LOW_Res data upload ####
Low_Res_IspH_AIIEIVDALDR_Spectra <- read.csv("F:/Projects/Proteomics/Zymomona/FAIMS/DataAnalysis/LibraryComparison/LowRes_LibrarySpectra_AIEIVDQALDR.csv", 
                                             header = TRUE, sep = ",", stringsAsFactors = FALSE)
Low_Res_IspH_AIIEIVDALDR_Spectra.RelaiveAbundance <- Low_Res_IspH_AIIEIVDALDR_Spectra
Low_Res_IspH_AIIEIVDALDR_Spectra.RelaiveAbundance$Relative.Abundance <- Low_Res_IspH_AIIEIVDALDR_Spectra.RelaiveAbundance$Intensity/(max(Low_Res_IspH_AIIEIVDALDR_Spectra.RelaiveAbundance$Intensity))
colnames(Low_Res_IspH_AIIEIVDALDR_Spectra.RelaiveAbundance) <- c("m.z", "Intensity", "Relative.Abundance")


Low_Res_IspG_ETDIGVTGGGGQGK_Spectra <- read.csv("F:/Projects/Proteomics/Zymomona/FAIMS/DataAnalysis/LibraryComparison/LowRes_LibrarySpectra_ETDIGVTGGGGQGK.csv", 
                                                header = TRUE, sep = ",", stringsAsFactors = FALSE)
Low_Res_IspG_ETDIGVTGGGGQGK_Spectra.RelativeAbundance <- Low_Res_IspG_ETDIGVTGGGGQGK_Spectra
Low_Res_IspG_ETDIGVTGGGGQGK_Spectra.RelativeAbundance$Relative.Abundance <- Low_Res_IspG_ETDIGVTGGGGQGK_Spectra.RelativeAbundance$Intensity/(max(Low_Res_IspG_ETDIGVTGGGGQGK_Spectra.RelativeAbundance$Intensity))
colnames(Low_Res_IspG_ETDIGVTGGGGQGK_Spectra.RelativeAbundance) <- c("m.z", "Intensity", "Relative.Abundance")


ggplot(PROSIT_IspH_AIIEIVDALDR_Spectra, aes(x=m.z, y= Relative.Abundance)) +
  #geom_line()+
  geom_linerange(aes(x=m.z, ymax=Relative.Abundance, ymin=0),
                 position = position_jitter(height = 0L, seed = 1L))+
  xlim(150,1250)+
  theme_classic()+
  labs(title = "PROSIT AIEIVDQDALDR Spectra", y = "Relative Abundance", x ="m/z")

ggplot(PROSIT_IspG_ETDIGVTGGGGQGK_Spectra, aes(x=m.z, y= Relative.Abundance)) +
  #geom_line()+
  geom_linerange(aes(x=m.z, ymax=Relative.Abundance, ymin=0),
                 position = position_jitter(height = 0L, seed = 1L))+
  xlim(150,1250)+
  theme_classic()+
  labs(title = "PROSIT ETDIGVTGGGGQGK Spectra", y = "Relative Abundance", x ="m/z")


######### Plotting Spectra ############

peptide_Library_fragment <- function(df_library){
  
  fragment <- numeric()
  for (i in 1:nrow(df_library)){
    if(df_library[i,2] != 0 ){
      fragment <- append(fragment,as.numeric(df_library[i,1]))
    }
  }
  
  return <- fragment
}


# df_library <- ProteinProspector_AIEIVQALDR
# df_prosit <- PROSIT_IspH_AIIEIVDALDR_Spectra
top10_ions <- function(df_library,df_prosit){
  
  # extract peptide fragments from protein prospect library
  fragment <- peptide_Library_fragment(df_library)
  
  # Print the name of the dataframe that will be compared to library peptides
  print(deparse(substitute(df)))
  
  # Create an empty vector to put the transitions that were found in the experimental spectra
  transitions <- numeric()
  
  # Ion comparison has to be exact, there for I will truncate m/z values before the decimal
  # and then turn them into characters for proper comparison
  df.mass.to.charge.char <- as.character(trunc(df_prosit$m.z))
  
  # iterate through library ions - fragment vector
  for(i in 1:length(fragment)){
    
    # truncate and turn fragment ion into a character
    frag_char <- as.character(trunc(fragment[i]))
    
    # assign a new variable for exact pattern matching
    pat <- paste0("^",frag_char,"$")
    
    # append matched ions to transition vector
    # transitions object will be populated with the index of where the match was found
    transitions <- append(transitions,grep(pat, 
                                           df.mass.to.charge.char))
  }
  
  mass.to.charge <- df_prosit$m.z[transitions]
  print(mass.to.charge)
  
  df_sort <-df_prosit[transitions,]
  df_sort <- df_sort[order(df_sort$Relative.Abundance, decreasing = TRUE),]
  
  df_top10 <- df_sort[1:10,]
  
  return <- df_top10
}

AIEIVDQALDR_top10 <- top10_ions(ProteinProspector_AIEIVQALDR, PROSIT_IspH_AIIEIVDALDR_Spectra)
ETDIGVTGGGQGK_top10 <- top10_ions(ProteinProspector_ETDIGVTGGGQGK, PROSIT_IspG_ETDIGVTGGGGQGK_Spectra)

final_list_using_Intensity_function <- function(df,transitions_vector){
  
  mass.to.charge <- df$m.z[transitions]
  print(mass.to.charge)
  
  # reduce to only the unique m/z
  unique.m.t.c <- unique(trunc(mass.to.charge))
  
  
  final_list <- vector()
  
  # The following will be used to match the m/z based on intensity. If there is more than one m/z
  # value being compared, the m/z with the greatest intensity is chosen. 
  for(i in 1:length(unique.m.t.c)){
    comparison_intensity <- which(trunc(mass.to.charge) %in% unique.m.t.c[i])
    
    if(length(comparison_intensity) <= 1){
      keep <- mass.to.charge[comparison_intensity]
      
      final_list <- append(final_list,keep)
      print(final_list)
      
      
    }else if(length(comparison_intensity) > 1){
      
      keep.intensity <- max(df$Intensity[which(df$m.z %in% mass.to.charge[comparison_intensity])])
      
      keep <- df$m.z[which(df$Intensity == keep.intensity)]
      
      final_list <- append(final_list, keep)
      print(final_list)
      
    }
    
  }
  
  return <- final_list
}

df_1 <- df_transitions
transitions_vector <- transitions
fragment_vector <- fragment

observe_ppm_forIons <- function(df_1, transitions_vector, fragment_vector) {
  mass.to.charge <- df_1$m.z
  print(mass.to.charge)
  
  # reduce to only the unique m/z
  unique.m.t.c <- unique(trunc(mass.to.charge))
  
  
  compiled_ppm_list <- vector(mode = "list", length = length(unique.m.t.c))
  # The following will be used to match the m/z based on intensity. If there is more than one m/z
  # value being compared, the m/z within the 10 ppm tolerance is chosen. 
  for(i in 1:length(unique.m.t.c)){
    
    comparison_intensity <- which(trunc(mass.to.charge) %in% unique.m.t.c[i])
    
    compiled_ppm <- vector()
    
    
    for(j in 1:length(comparison_intensity)){
      theoretical <- fragment_vector[which(trunc(fragment_vector) %in% unique.m.t.c[i])]
      observed <- mass.to.charge[comparison_intensity[j]]
      
      ppm_error <- abs(observed-theoretical)/theoretical*1e6
      
      compiled_ppm <- append(compiled_ppm, ppm_error)
      
      
    }
    
    ppm_Df <- data.frame(matrix(NA, nrow = length(comparison_intensity), ncol = 4))
    names(ppm_Df) <- c("m.z", "Intensity", "Relative.Abundance", "ppm")
    
    
    ppm_Df$m.z <- mass.to.charge[comparison_intensity]
    ppm_Df$Intensity <- df_1$Intensity[which(df_1$m.z %in% mass.to.charge[comparison_intensity])]
    ppm_Df$Relative.Abundance <- df_1$Relative.Abundance[which(df_1$m.z %in% mass.to.charge[comparison_intensity])]
    ppm_Df$ppm <- compiled_ppm
    
    compiled_ppm_list[[i]] <- ppm_Df
   
  }
  
    names(compiled_ppm_list) <- unique.m.t.c
    
    return(compiled_ppm_list)

  }


df_1 <- df_transitions
transitions_vector <- transitions
fragment_vector <- fragment

final_list_using_infusion_function <- function(df_1, transitions_vector, fragment_vector){
  
  mass.to.charge <- df_1$m.z
  print(mass.to.charge)
  
  # reduce to only the unique m/z
  unique.m.t.c <- unique(trunc(mass.to.charge))
  
  list_to_chooseFrom <- observe_ppm_forIons(df_1, transitions_vector, fragment_vector)
  
  compiled_peptides <- vector()
   
   for (i in 1:length(list_to_chooseFrom)) {
     extract_df <- list_to_chooseFrom[[i]]
     
     extract_df_lessthan15ppm <- extract_df[which(extract_df$ppm < 15),]
     
     extract_df_lessthan15ppm_highAbundance <- extract_df_lessthan15ppm[which(extract_df_lessthan15ppm$Relative.Abundance == 
                                                                                max(extract_df_lessthan15ppm$Relative.Abundance)),]
     
     compiled_peptides <- append(compiled_peptides, extract_df_lessthan15ppm_highAbundance$m.z)
     
     }
     
    return(compiled_peptides)
   }


df_1 <- df_transitions
transitions_vector <- transitions
fragment_vector <- fragment
 
# match fragments for high resolution data
final_list_using_ppm_function <- function(df_1,transitions_vector,fragment_vector){
  
  mass.to.charge <- df_1$m.z
  print(mass.to.charge)
  
  # reduce to only the unique m/z
  unique.m.t.c <- unique(trunc(mass.to.charge))
  
  
  final_list <- vector()
  
  # The following will be used to match the m/z based on intensity. If there is more than one m/z
  # value being compared, the m/z within the 10 ppm tolerance is chosen. 
  for(i in 1:length(unique.m.t.c)){
    
    comparison_intensity <- which(trunc(mass.to.charge) %in% unique.m.t.c[i])
       
      within_ppm <- vector()
      
      for(j in 1:length(comparison_intensity)){
        
        theoretical <- fragment_vector[which(trunc(fragment_vector) %in% unique.m.t.c[i])]
        observed <- mass.to.charge[comparison_intensity[j]]
        
        ppm_error <- abs(observed-theoretical)/theoretical*1e6
        
        if(ppm_error <= 10){
          
          within_ppm <- append(within_ppm, ppm_error)
           keep <- observed
           
           final_list <- append(final_list, keep)
           
           print(final_list)
        }
      }
  }
  
  df_final_list <- df_transitions[which(df_transitions$m.z %in% final_list),]
    return <- final_list
}

# match fragments for low resolution data
final_list_using_DAwindow_function <- function(df_1,transitions_vector,fragment_vector){
  mass.to.charge <- df_1$m.z[transitions_vector]
  #print(mass.to.charge)
  
  # reduce to only the unique m/z
  unique.m.t.c <- unique(trunc(mass.to.charge))
  
  
  final_list <- vector()
  
  # The following will be used to match the m/z based on intensity. If there is more than one m/z
  # value being compared, the m/z with the greatest intensity is chosen. 
  for(i in 1:length(unique.m.t.c)){
    
    comparison_intensity <- which(trunc(mass.to.charge) %in% unique.m.t.c[i])
      
      #within_Da_window <- vector()
      
      for(j in 1:length(comparison_intensity)){
        
       theoretical = fragment_vector[which(trunc(fragment_vector) %in% unique.m.t.c[i])]
       theoretical_lowEnd = theoretical-0.35
       theoretical_highEnd = theoretical+0.35
       
       observed = mass.to.charge[comparison_intensity[j]]
        
       if(between(observed,theoretical_lowEnd,theoretical_highEnd))
       {
         #within_Da_window <- append(within_Da_window, observed)
         
          keep <- observed
          
          final_list <- append(final_list, keep)
          
          print(final_list)
          
          #rm(observed)
       #}
      }
    }
  }
  return <- final_list
}

 # df_library = AIEIVDQALDR_top10
 # df = FAIMS_AIIEIVDALDR_Spectra_RelativeAbundance
 df_library = ETDIGVTGGGQGK_top10
 df = noFAIMS_ETDIGVTGGGQGK_Spectra_RelativeAbundance
spectra_plotting <- function(df_library, df){
  # extract the library peptide ions into object fragment
  fragment <- sort(peptide_Library_fragment(df_library))
  
  
  # Print the name of the dataframe that will be compared to library peptides
  print(deparse(substitute(df)))
  # Create an empty vector to put the transitions that were found in the experimental spectra
  transitions <- numeric()
  
  # Ion comparison has to be exact, there for I will truncate m/z values before the decimal
  # and then turn them into characters for proper comparison
  df.mass.to.charge.char <- as.character(trunc(df$m.z))
  
  # iterate through library ions - fragment vector
  for(i in 1:length(fragment)){
    
    # truncate and turn fragment ion into a character
    frag_char <- as.character(trunc(fragment[i]))
    
    # assign a new variable for exact pattern matching
    pat <- paste0("^",frag_char,"$")
    
    # append matched ions to transition vector
    # transitions object will be populated with the index of where the match was found
    transitions <- append(transitions,grep(pat, 
                                           df.mass.to.charge.char))
    
    df_transitions <- df[transitions,]
  }
  
  if(length(grep("High",deparse(substitute(df)))) != 0){
    print("use 10 ppm window")
    final_list <- final_list_using_ppm_function(df, transitions, fragment)
  }else if(length(grep("Low",deparse(substitute(df)))) != 0){
    print("use 0.35 Da window")
    final_list <- final_list_using_DAwindow_function(df, transitions, fragment)
  }else if(length(grep("PROSIT",deparse(substitute(df)))) != 0){
    final_list <- df$m.z[transitions]
  }else if(length(grep("FAIMS",deparse(substitute(df)))) != 0){
    print("use ppm and intensitywindow")
    final_list <- final_list_using_infusion_function(df_transitions, transitions, fragment)
  }
  
  
  # Differentiate between the ions that are associatted to the peptide and those that are background.
  # Color the peptide ions red and the background black
  ion_type <- vector()

  
  for (i in 1:nrow(df)) {
    if(length(which(i %in% which(df$m.z %in% final_list))) == 1){
      ion_type <- append(ion_type, "peptide")
    }else{
      ion_type <- append(ion_type, "background")
    }
  }
  
  background <- df[grep("background",ion_type),]
  peptide <- df[grep("peptide",ion_type),]
  
  # produce the ggplot that will be returned with the function
  spectra <- ggplot(background, aes(x=m.z, y= Relative.Abundance)) +
    #geom_line()+
    geom_linerange(aes(x=m.z, ymax=Relative.Abundance, ymin=0),
                   position = position_jitter(height = 0L, seed = 1L))+
    geom_linerange(data = peptide, 
                   aes(x=m.z, ymax = Relative.Abundance, ymin =0), 
                   position = position_jitter(height = 0L, seed = 1L), color = "red")+
    xlim(150,1250)+
    theme_classic()
  
  return <- spectra
  
}

ggplot(df_transitions, aes(x=m.z, y= Intensity))+
  geom_linerange(aes(x=m.z, ymax=Intensity, ymin=0),position = position_jitter(height = 0L, seed = 1L))

##### Plot spectra #####
prosit_AIEIVDQALDR <- spectra_plotting(AIEIVDQALDR_top10, PROSIT_IspH_AIIEIVDALDR_Spectra)
PA <- prosit_AIEIVDQALDR + labs(title = "Prosit Theoretical Spectra - AIEVDQALDR", y = "Relative Abundance", x = "m/z")

high_res_AIEIDVDQDALDR <- spectra_plotting(AIEIVDQALDR_top10, High_Res_IspH_AIIEIVDALDR_Spectra.RelativeAbundance)
HA <- high_res_AIEIDVDQDALDR + labs(title = "High Resolution AIEIVDQDALDR Spectra", y = "Relative Abundance", x ="m/z")

low_res_AIEIDVDQDALDR <- spectra_plotting(AIEIVDQALDR_top10, Low_Res_IspH_AIIEIVDALDR_Spectra.RelaiveAbundance)
LA <- low_res_AIEIDVDQDALDR + labs(title = "Low Resolution AIEIVDQDALDR Spectra", y = "Relative Abundance", x ="m/z")
#labs(title = "Low Resolution ETDIGVTGGGGQGK Spectra", y = "Relative Abundance", x ="m/z")

prosit_ETDIGVTGGGGQGK <- spectra_plotting(ETDIGVTGGGQGK_top10, PROSIT_IspG_ETDIGVTGGGGQGK_Spectra)
PE <- prosit_ETDIGVTGGGGQGK + labs(title = "Prosit Theoretical Spectra - ETDIGVTGGGQGK", y= "Relative Abundance", x="m/z")

high_res_ETDIGVTGGGGQGK <- spectra_plotting(ETDIGVTGGGQGK_top10, High_Res_IspG_ETDIGVTGGGGQGK_Spectra.RelativeAbundance)
HE <- high_res_ETDIGVTGGGGQGK + labs(title = "High Resolution ETDIGVTGGGGQGK Spectra", y = "Relative Abundance", x ="m/z")

low_res_ETDIGVTGGGGQGK <- spectra_plotting(ETDIGVTGGGQGK_top10, Low_Res_IspG_ETDIGVTGGGGQGK_Spectra.RelativeAbundance)
LE <- low_res_ETDIGVTGGGGQGK + labs(title = "Low Resolution ETDIGVTGGGGQGK Spectra", y = "Relative Abundance", x ="m/z")


##### Infusion Comparison to Library ######
# The following functin infusion_spectra_prep is meant to format the infusion spectra for generating 
# the spectral plots

Infusion_df <- FAIMS_ETDIGVTGGGQGK_Spectra
max_Intensity <- 200000

infusion_spectra_prep <- function(Infusion_df, max_Intensity){
  colnames(Infusion_df) <- c("m.z", "Intensity")
  
  row_sub <- apply(Infusion_df, 1, function(row) all(row !=0))
  Infusion_df <- Infusion_df[row_sub,]
  
  Infusion_df <- Infusion_df[order(Infusion_df$m.z),]
  
  Infusion_df_RelativeAbundance <- Infusion_df[(Infusion_df$Intensity < max_Intensity),]
  
  Infusion_df_RelativeAbundance$Relative.Abundance <- Infusion_df_RelativeAbundance$Intensity/(max(Infusion_df_RelativeAbundance$Intensity))
  
  return(Infusion_df_RelativeAbundance)
}


##### Infusion data upload #####
#AIEIVDQALDR
#FAIMS_AIEIVDALDR_old_v = "F:/Projects/Proteomics/Zymomona/FAIMS/DataAnalysis/LibraryComparison/20200126_ZymoFAIMS_IW2_R500K_MI502_AGC1e06_2.csv"
FAIMS_AIIEIVDALDR_Spectra <- read.csv("P:/EAT_20190926_Zymomona/ParameterComparisonSpectra/20200126_ZymoFAIMS_IW_0_6_R240K_MI502_AGC1e06_1_AIEIVDQALDR_R.csv", 
                                      header = TRUE, sep = ",", stringsAsFactors = FALSE)

FAIMS_AIIEIVDALDR_Spectra_RelativeAbundance <- infusion_spectra_prep(FAIMS_AIIEIVDALDR_Spectra, 200000)

FAIMS_AIEIDVDQDALDR <- spectra_plotting(AIEIVDQALDR_top10, FAIMS_AIIEIVDALDR_Spectra_RelativeAbundance)
A1 <- FAIMS_AIEIDVDQDALDR + ylim(0,1) + labs(title = "FAIMS AIEIVDQDALDR Spectra", y = "Relative Abundance", x ="m/z")

#with out FAIMS
noFAIMS_AIIEIVDALDR_Spectra <- read.csv("P:/EAT_20190926_Zymomona/ParameterComparisonSpectra/20200125_Zymo_IW_0_6_R240K_MI502_AGC1e06_1_AIEIVDQALDR_R.csv", 
                                        header = TRUE, sep = ",", stringsAsFactors = FALSE)

noFAIMS_AIIEIVDALDR_Spectra_RelativeAbundance <- infusion_spectra_prep(noFAIMS_AIIEIVDALDR_Spectra, 200000)

noFAIMS_AIEIDVDQDALDR <- spectra_plotting(AIEIVDQALDR_top10, noFAIMS_AIIEIVDALDR_Spectra_RelativeAbundance)
A2 <- noFAIMS_AIEIDVDQDALDR + ylim(0,1) + labs(title = "no FAIMS AIEIVDQDALDR Spectra", y = "Relative Abundance", x ="m/z")

# ETDIGVTGGGQGK #

#FAIMS
#FAIMS_ETDIGVTGGGQGK_old_v = "F:/Projects/Proteomics/Zymomona/FAIMS/DataAnalysis/LibraryComparison/20200126_ZymoFAIMS_IW2_R500K_MI502_AGC1e06_2_ETDIGVTGGGQGK.csv"
FAIMS_ETDIGVTGGGQGK_Spectra <- read.csv("P:/EAT_20190926_Zymomona/ParameterComparisonSpectra/20200126_ZymoFAIMS_IW_0_6_R240K_MI502_AGC1e06_1_ETDIGVTGGGQGK_R.csv", 
                                        header = TRUE, sep = ",", stringsAsFactors = FALSE)

FAIMS_ETDIGVTGGGQGK_Spectra_RelativeAbundance <- infusion_spectra_prep(FAIMS_ETDIGVTGGGQGK_Spectra, 200000)

FAIMS_ETDIGVTGGGQGK <- spectra_plotting(ETDIGVTGGGQGK_top10, FAIMS_ETDIGVTGGGQGK_Spectra_RelativeAbundance)
E1 <- FAIMS_ETDIGVTGGGQGK + ylim(0,1) + labs(title = "FAIMS ETDIGVTGGGQGK Spectra", y = "Relative Abundance", x ="m/z")

#with out FAIMS
noFAIMS_ETDIGVTGGGQGK_Spectra <- read.csv("P:/EAT_20190926_Zymomona/ParameterComparisonSpectra/20200125_Zymo_IW_0_6_R240K_MI502_AGC1e06_2_ETDIGVTGGGQGK_R.csv", 
                                          header = TRUE, sep = ",", stringsAsFactors = FALSE)

noFAIMS_ETDIGVTGGGQGK_Spectra_RelativeAbundance <- infusion_spectra_prep(noFAIMS_ETDIGVTGGGQGK_Spectra, 200000)

noFAIMS_ETDIGVTGGGGQGK <- spectra_plotting(ETDIGVTGGGQGK_top10, noFAIMS_ETDIGVTGGGQGK_Spectra_RelativeAbundance)
E2 <- noFAIMS_ETDIGVTGGGGQGK + ylim(0,1) + labs(title = "no FAIMS ETDIGVTGGGQGK Spectra", y = "Relative Abundance", x ="m/z")

# Use the ggplotGrob function to prepare plots to be arranged for pdf printing
PAgrob <- ggplotGrob(PA)
HAgrob <- ggplotGrob(HA)
LAgrob <- ggplotGrob(LA)

PEgrob <- ggplotGrob(PE)
HEgrob <- ggplotGrob(HE)
LEgrob <- ggplotGrob(LE)

A1grob <- ggplotGrob(A1)
A2grob <- ggplotGrob(A2)
E1grob <- ggplotGrob(E1)
E2grob <- ggplotGrob(E2)


grid.arrange(PAgrob, PEgrob, HAgrob, HEgrob, LAgrob, LEgrob, ncol =2)

grid.arrange(A1grob, E1grob, A2grob, E2grob, ncol = 2)


pdf("F:/Projects/Proteomics/Zymomona/FAIMS/Figures/FromR/SpectraComparisons_PROSIT_HIGHRES_LOWRES_IW6.pdf")
grid.arrange(PAgrob, PEgrob, HAgrob, HEgrob, LAgrob, LEgrob, ncol =2)
dev.off()


ggsave("F:/Projects/Proteomics/Zymomona/FAIMS/Figures/FromR/SpectraComparisons_PROSIT_HIGHRES_LOWRES_IW6.pdf",
       grid.arrange(PAgrob, PEgrob, HAgrob, HEgrob, LAgrob, LEgrob, ncol =2),
       width = 10,
       height = 11.5,
       units = "in")

ggsave("F:/Projects/Proteomics/Zymomona/FAIMS/Figures/FromR/SpectraComparisons_FAIMS_noFAIMS_IW6.pdf",
       grid.arrange(A1grob, E1grob, A2grob, E2grob, ncol =2),
       width = 10,
       height = 8.5,
       units = "in")

ggsave("F:/Projects/Proteomics/Zymomona/FAIMS/Figures/FromR/SpectraComparisons_all_top10_IW6.pdf",
       grid.arrange(PAgrob, PEgrob, HAgrob, HEgrob, LAgrob, LEgrob, A1grob, E1grob, A2grob, E2grob, ncol =2),
       width = 10,
       height = 14.5,
       units = "in")

#### TIC Explained ######
# To explain the TIC associated to the peptide fragments the function tic_explained_function

 df_proteinProspect <- ETDIGVTGGGQGK_top10
 df <- FAIMS_ETDIGVTGGGQGK_Spectra_RelativeAbundance

tic_explained_function_LC <- function(df_proteinProspect, df){
  
  #fragment <- df_proteinProspect$m.z
  
  transitions <- numeric()
  
  df.mass.to.charge.char <- as.character(sort(trunc(df_proteinProspect$m.z)))
  
  
  for(i in 1:length(df.mass.to.charge.char)){
    
    frag_char <- df.mass.to.charge.char[i]
    
    pat <- paste0("^",frag_char,"$")
    
    transitions <- append(transitions,grep(pat,as.character(trunc(df$m.z))))
  }
  
  mass.to.charge <- df$m.z[transitions]
  print(mass.to.charge)
  
  unique.m.t.c <- unique(trunc(mass.to.charge))
  
  final_list <- vector()

  for(i in 1:length(unique.m.t.c)){
    comparison_intensity <- which(trunc(mass.to.charge) %in% unique.m.t.c[i])

    if(length(comparison_intensity) <= 1){
      keep <- mass.to.charge[comparison_intensity]

      final_list <- append(final_list,keep)
      print(final_list)


    }else if(length(comparison_intensity) > 1){

      keep.intensity <- max(df$Intensity[which(df$m.z %in% mass.to.charge[comparison_intensity])])

      keep <- df$m.z[which(df$Intensity == keep.intensity)]

      final_list <- append(final_list, keep)
      print(final_list)
    }

  }
  
  intensity_peptide <- df$Intensity[which(df$m.z %in% final_list)]
  intensity_peptide_sum <- sum(intensity_peptide)
  intensity_TIC <- sum(df$Intensity)
  
  Tic_explained <- intensity_peptide_sum/intensity_TIC*100
  
  return(Tic_explained)
  
}

tic_explained_function_infusion <- function(df_proteinProspect, df){
  
  #fragment <- df_proteinProspect$m.z
  
  transitions <- numeric()
  
  df.mass.to.charge.char <- as.character(sort(trunc(df_proteinProspect$m.z)))
  
  
  for(i in 1:length(df.mass.to.charge.char)){
    
    frag_char <- df.mass.to.charge.char[i]
    
    pat <- paste0("^",frag_char,"$")
    
    transitions <- append(transitions,grep(pat,as.character(trunc(df$m.z))))
  }
  
  mass.to.charge <- df$m.z[transitions]
  print(mass.to.charge)
  
  unique.m.t.c <- unique(trunc(mass.to.charge))
  
  final_list <- vector()
  contaminant <- vector()
  
  for(i in 1:length(unique.m.t.c)){
    comparison_intensity <- which(trunc(mass.to.charge) %in% unique.m.t.c[i])
    
    if(length(comparison_intensity) <= 1){
      keep <- mass.to.charge[comparison_intensity]
      
      final_list <- append(final_list,keep)
      print(final_list)
      
      
    }else if(length(comparison_intensity) > 1){
      
      keep.intensity <- max(df$Intensity[which(df$m.z %in% mass.to.charge[comparison_intensity])])
      
      keep <- df$m.z[which(df$Intensity == keep.intensity)]
      remove_contaminant.location <- which(df$Intensity == keep.intensity)+1
      
      remove_contaminant <- df$m.z[remove_contaminant.location]
      
      final_list <- append(final_list, keep)
      contaminant <- append(contaminant, remove_contaminant)
      print(final_list)
    }
    
  }
  
  intensity_cleaned <- df$Intensity[-which(df$m.z %in% contaminant)]
  intensity_peptide <- df$Intensity[which(df$m.z %in% final_list)]
  intensity_peptide_sum <- sum(intensity_peptide)
  intensity_TIC <- sum(intensity_cleaned)
  
  Tic_explained <- intensity_peptide_sum/intensity_TIC*100
  
  return(Tic_explained)
  
}

Tic_explained_High_Res_AIEIDVQALDR <- tic_explained_function(AIEIVDQALDR_top10, High_Res_IspH_AIIEIVDALDR_Spectra)

Tic_explained_Low_Res_AIEIDVQALDR <- tic_explained_function(AIEIVDQALDR_top10,Low_Res_IspH_AIIEIVDALDR_Spectra)

Tic_explained_High_Res_ETDIGVTGGGQGK <- tic_explained_function(ETDIGVTGGGQGK_top10, High_Res_IspG_ETDIGVTGGGGQGK_Spectra)

Tic_explained_Low_Res_ETDIGVTGGGQGK <- tic_explained_function(ETDIGVTGGGQGK_top10, Low_Res_IspG_ETDIGVTGGGGQGK_Spectra)

### AIEIVDQALDR ###
colnames(FAIMS_AIIEIVDALDR_Spectra_RelativeAbundance) <- c("m.z", "Intensity", "Relative.Abundance")
Tic_explained_FAIMS_AIEIDVQALDR <- tic_explained_function(AIEIVDQALDR_top10,FAIMS_AIIEIVDALDR_Spectra_RelativeAbundance)

colnames(noFAIMS_AIIEIVDALDR_Spectra_RelativeAbundance) <- c("m.z", "Intensity", "Relative.Abundance")
Tic_explained_noFAIMS_AIEIDVQALDR <- tic_explained_function(AIEIVDQALDR_top10, noFAIMS_AIIEIVDALDR_Spectra_RelativeAbundance)

### ETDIGVTGGGQGK ###
colnames(FAIMS_ETDIGVTGGGQGK_Spectra_RelativeAbundance) <- c("m.z", "Intensity", "Relative.Abundance")
Tic_explained_FAIMS_ETDIGVTGGGQGK <- tic_explained_function(ETDIGVTGGGQGK_top10,FAIMS_ETDIGVTGGGQGK_Spectra_RelativeAbundance)

colnames(noFAIMS_ETDIGVTGGGQGK_Spectra_RelativeAbundance) <- c("m.z", "Intensity", "Relative.Abundance")
Tic_explained_noFAIMS_ETDIGVTGGGQGK <- tic_explained_function(ETDIGVTGGGQGK_top10, noFAIMS_ETDIGVTGGGQGK_Spectra_RelativeAbundance)



##### Dot Product ####
 df_top10_comparison <- AIEIVDQALDR_top10
 df <- PROSIT_IspH_AIIEIVDALDR_Spectra

peptide_ions_function <- function(df_top10_comparison, df, df_name){
  
  print(df_name)
  
  transitions <- numeric()
  
  df.mass.to.charge.char <- as.character(sort(trunc(df_top10_comparison$m.z)))
  
  
  for(i in 1:length(df.mass.to.charge.char)){
    
    frag_char <- df.mass.to.charge.char[i]
    
    pat <- paste0("^",frag_char,"$")
    
    transitions <- append(transitions,grep(pat,as.character(trunc(df$m.z))))
  }
  
  
  
  #mass.to.charge <- df$m.z[transitions]
  #print(mass.to.charge)
  
  fragment <- sort(df_top10_comparison$m.z)
  
  if(length(grep("High",df_name) != 0)){
    
    print("use 10 ppm window")
      final_list <- final_list_using_ppm_function(df, transitions,fragment)
      peptide <- df[which(df$m.z %in% final_list),]
      
  }
  
  if(length(grep("Low",df_name) != 0)){
    
    print("use 0.35 Da window")
    final_list <- final_list_using_DAwindow_function(df, transitions, fragment)
    peptide <- df[which(df$m.z %in% final_list),]
    
  }
  
  if(length(grep("FAIMS",df_name)) != 0){
    
    print("use 10 ppm window and recalculate relative abundance")
    final_list <- final_list_using_ppm_function(df, transitions,fragment)
    
    unique.m.t.c <- unique(trunc(final_list))
    
    
    reduced_list <- vector()
    
    # The following will be used to match the m/z based on intensity. If there is more than one m/z
    # value being compared, the m/z with the greatest intensity is chosen. 
    for(i in 1:length(unique.m.t.c)){
      comparison_intensity <- which(trunc(final_list) %in% unique.m.t.c[i])
      
      if(length(comparison_intensity) <= 1){
        keep <- final_list[comparison_intensity]
        
        reduced_list <- append(reduced_list,keep)
        print(reduced_list)
        
        
      }else if(length(comparison_intensity) > 1){
        
        
        keep.intensity <- max(df$Intensity[which(df$m.z %in% final_list[comparison_intensity])])
        
        keep <- df$m.z[which(df$Intensity == keep.intensity)]
        
        reduced_list <- append(reduced_list, keep)
        print(reduced_list)
        
      }
      
    }
    
    peptide <- df[which(df$m.z %in% reduced_list),]
    
    max_intensity <- max(peptide$Intensity)
    
    peptide$Recalculated.Relative.Abundance <- peptide$Intensity/max_intensity
    
  }
  if(length(grep("PROSIT",df_name)) != 0){
    peptide <- df[which(df$m.z %in% fragment),]
  }
  

  
  
  return(peptide)
  
}



df_proteinProspect <- AIEIVDQALDR_top10
df <- noFAIMS_AIIEIVDALDR_Spectra_RelativeAbundance

intensity_function <- function(df_proteinProspect, df){
  
  transitions <- numeric()
  
  df.mass.to.charge.char <- as.character(sort(trunc(df_proteinProspect$m.z)))
  
  
  for(i in 1:length(df.mass.to.charge.char)){
    
    frag_char <- df.mass.to.charge.char[i]
    
    pat <- paste0("^",frag_char,"$")
    
    transitions <- append(transitions,grep(pat,as.character(trunc(df$m.z))))
  }
  
  
  mass.to.charge <- df$m.z[transitions]
  #print(mass.to.charge)
  
  unique.m.t.c <- unique(trunc(mass.to.charge))
  
  
  final_list <- vector()
  
  for(i in 1:length(unique.m.t.c)){
    comparison_intensity <- which(trunc(mass.to.charge) %in% unique.m.t.c[i])
    
    if(length(comparison_intensity) <= 1){
      keep <- mass.to.charge[comparison_intensity]
      
      final_list <- append(final_list,keep)
      #print(final_list)
      
      
    }else if(length(comparison_intensity) > 1){
      
      keep.intensity <- max(df$Relative.Abundance[which(df$m.z %in% mass.to.charge[comparison_intensity])])
      
      keep <- df$m.z[which(df$Relative.Abundance == keep.intensity)]
      
      final_list <- append(final_list, keep)
      #print(final_list)
    }
    
  }
  
  intensity_peptide <- df[which(df$m.z %in% final_list),ncol(df)]
  #intensity_peptide_sum <- sum(intensity_peptide)
  
  return(intensity_peptide)
  
}

  df_1 <- PROSIT_IspH_AIIEIVDALDR_Spectra
  df_2 <- noFAIMS_AIIEIVDALDR_Spectra_RelativeAbundance
  df_top10 <- AIEIVDQALDR_top10

DP_function <- function(df_1, df_2, df_top10, name1, name2){
  
  df1_peptide <- peptide_ions_function(df_top10,df_1,name1) 
  
  df2_peptide <- peptide_ions_function(df_top10, df_2,name2)
  
  if(nrow(df1_peptide) >= nrow(df2_peptide)){
    spec1 = df1_peptide
    spec2 = df2_peptide
  }else{
    spec1 = df2_peptide
    spec2 = df1_peptide
  }
  
  #sum_intensity_combined = vector()
  sum_Relative.Abundance_combined = vector()
  not_included = vector()
  
  for(i in 1:nrow(spec1)){
    
    print(spec1$m.z[i])
    
    if(is.na(trunc(spec2$m.z[i]))){
      print(paste0("i = ",i," Ion NOt here"))
      print(spec1$m.z[i])
      print(spec2$m.z[i])
    }
    
    if(length(which(trunc(spec1$m.z[i]) %in% trunc(spec2$m.z))) != 0){
      
      #spec1_intensity = spec1$Intensity[i]
      spec1_Relative.Abundance = spec1[,ncol(spec1)][i]
      
      #spec2_intensity = spec2$Intensity[which(trunc(spec2$m.z) == trunc(spec1$m.z[i]))]
      spec2_Relative.Abundance = spec2[,ncol(spec2)][which(trunc(spec2$m.z) == trunc(spec1$m.z[i]))]
      
      
      #combined_intensity = spec1_intensity*spec2_intensity
      combined_Relative.Abundance = spec1_Relative.Abundance*spec2_Relative.Abundance
      
      # sum_intensity_combined = append(sum_intensity_combined, combined_intensity)
      # print(sum_intensity_combined)
      
      sum_Relative.Abundance_combined = append(sum_Relative.Abundance_combined, combined_Relative.Abundance)
      print(sum_Relative.Abundance_combined)
      
      
    } else{
      
      print(paste0("i = ",i," Ion NOt here"))
      
      #spec1_intensity = spec1$Intensity[i]
      spec1_Relative.Abundance = spec1[,ncol(spec1)][i]
      
      #spec2_intensity = spec2$Intensity[which(trunc(spec2$m.z) == trunc(spec1$m.z[i]))]
      spec2_Relative.Abundance = 0
      
      #combined_intensity = spec1_intensity*spec2_intensity
      combined_Relative.Abundance = spec1_Relative.Abundance*spec2_Relative.Abundance
      
      sum_Relative.Abundance_combined = append(sum_Relative.Abundance_combined, combined_Relative.Abundance)
      print(sum_Relative.Abundance_combined)
      
    }
  }
  
  
  I_spec1 <- intensity_function(df_top10, df1_peptide)
  I_spec2 <- intensity_function(df_top10, df2_peptide)
  
  numerator <- sum(sum_Relative.Abundance_combined)
  denomenator <- sum(I_spec1^2)*sum(I_spec2^2)
  
  DP <- numerator/sqrt(denomenator)
  
  return(DP)
  
}

DP_function(PROSIT_IspH_AIIEIVDALDR_Spectra,High_Res_IspH_AIIEIVDALDR_Spectra.RelativeAbundance, AIEIVDQALDR_top10, "PROSIT_IspH_AIIEIVDALDR_Spectra", "High_Res_IspH_AIIEIVDALDR_Spectra.RelativeAbundance")

DP_function(High_Res_IspH_AIIEIVDALDR_Spectra.RelativeAbundance,High_Res_IspH_AIIEIVDALDR_Spectra.RelativeAbundance, AIEIVDQALDR_top10)
DP_function(High_Res_IspH_AIIEIVDALDR_Spectra.RelativeAbundance,PROSIT_IspH_AIIEIVDALDR_Spectra, AIEIVDQALDR_top10,"High_Res_IspH_AIIEIVDALDR_Spectra.RelativeAbundance","PROSIT_IspH_AIIEIVDALDR_Spectra")
DP_function(High_Res_IspH_AIIEIVDALDR_Spectra.RelativeAbundance,Low_Res_IspH_AIIEIVDALDR_Spectra.RelaiveAbundance, AIEIVDQALDR_top10)
DP_function(High_Res_IspH_AIIEIVDALDR_Spectra.RelativeAbundance,FAIMS_AIIEIVDALDR_Spectra_RelativeAbundance, AIEIVDQALDR_top10)


DP_function(PROSIT_IspG_ETDIGVTGGGGQGK_Spectra,High_Res_IspG_ETDIGVTGGGGQGK_Spectra.RelativeAbundance, ETDIGVTGGGQGK_top10, "PROSIT_IspG_ETDIGVTGGGGQGK_Spectra", "High_Res_IspG_ETDIGVTGGGGQGK_Spectra.RelativeAbundance")

DP_function(High_Res_IspG_ETDIGVTGGGGQGK_Spectra.RelativeAbundance,noFAIMS_ETDIGVTGGGQGK_Spectra_RelativeAbundance, ETDIGVTGGGQGK_top10)


df.list <- list(PROSIT_IspH_AIIEIVDALDR_Spectra,High_Res_IspH_AIIEIVDALDR_Spectra.RelativeAbundance, Low_Res_IspH_AIIEIVDALDR_Spectra.RelaiveAbundance, 
                FAIMS_AIIEIVDALDR_Spectra_RelativeAbundance, noFAIMS_AIIEIVDALDR_Spectra_RelativeAbundance)
names(df.list) <- c("PROSIT_IspH_AIIEIVDALDR_Spectra", "High_Res_IspH_AIIEIVDALDR_Spectra.RelativeAbundance", "Low_Res_IspH_AIIEIVDALDR_Spectra.RelaiveAbundance",
                    "FAIMS_AIIEIVDALDR_Spectra_RelativeAbundance","noFAIMS_AIIEIVDALDR_Spectra_RelativeAbundance")

df.list_ETDIGVTGGGGQGK <- list(PROSIT_IspG_ETDIGVTGGGGQGK_Spectra,High_Res_IspG_ETDIGVTGGGGQGK_Spectra.RelativeAbundance, Low_Res_IspG_ETDIGVTGGGGQGK_Spectra.RelativeAbundance, 
                               FAIMS_ETDIGVTGGGQGK_Spectra_RelativeAbundance, noFAIMS_ETDIGVTGGGQGK_Spectra_RelativeAbundance)
names(df.list_ETDIGVTGGGGQGK) <- list("PROSIT_IspG_ETDIGVTGGGGQGK_Spectra","High_Res_IspG_ETDIGVTGGGGQGK_Spectra.RelativeAbundance", "Low_Res_IspG_ETDIGVTGGGGQGK_Spectra.RelativeAbundance", 
                               "FAIMS_ETDIGVTGGGQGK_Spectra_RelativeAbundance", "noFAIMS_ETDIGVTGGGQGK_Spectra_RelativeAbundance")


df_matrix <- matrix(nrow = 5, ncol = 5)
df_matrix_ETDIGVTGGGGQGK <- matrix(nrow = 5, ncol = 5)

for(i in 1:length(df.list)){
  df_n <- df.list[[i]]
  df_name <- names(df.list[i])
  for (j in 1:length(df.list)) {
    dot_product <- DP_function(df_n,df.list[[j]], AIEIVDQALDR_top10, df_name, names(df.list[j]))
    df_matrix[i,j] <- dot_product
  }
}

for(i in 1:length(df.list_ETDIGVTGGGGQGK)){
  df_n <- df.list_ETDIGVTGGGGQGK[[i]]
  df_name <- names(df.list_ETDIGVTGGGGQGK[i])
  for (j in 1:length(df.list_ETDIGVTGGGGQGK)) {
    dot_product = DP_function(df_n,df.list_ETDIGVTGGGGQGK[[j]], ETDIGVTGGGQGK_top10,df_name,names(df.list_ETDIGVTGGGGQGK[j]))
    df_matrix_ETDIGVTGGGGQGK[i,j] <- dot_product
  }
}

colnames(df_matrix) <- c("PROSIT","High Res", "Low Res", "Infusion FAIMS", "Infusion")
rownames(df_matrix) <- c("PROSIT","High Res", "Low Res", "Infusion FAIMS", "Infusion")

colnames(df_matrix_ETDIGVTGGGGQGK) <- c("PROSIT","High Res", "Low Res", "Infusion FAIMS", "Infusion")
rownames(df_matrix_ETDIGVTGGGGQGK) <- c("PROSIT","High Res", "Low Res", "Infusion FAIMS", "Infusion")

library(RColorBrewer)
#cols <- brewer.pal(4,"BrBG")
scaleRYG <- colorRampPalette(c("#F3A5BF","#15688E"), space = "rgb")(100)



pdf("F:/Projects/Proteomics/Zymomona/FAIMS/Figures/FromR/dotProduct_AIEIVDQALDR_ETDIGVTGGGQGK_IW6.pdf",height = 10, width = 15)
par(mfrow = c(1,2))
corrplot(df_matrix, type = "lower", method = "color", #col = scaleRYG, 
         addCoef.col = "white",
         tl.col = "black", tl.srt = 45,
         cl.lim = c(0,1))

corrplot(df_matrix_ETDIGVTGGGGQGK, type = "lower", method = "color", #col = scaleRYG, 
         addCoef.col = "white",
         tl.col = "black", tl.srt = 45,
         cl.lim = c(0,1))
dev.off()

#ggcorrplot(df_matrix, type = "lower", outline.color = "white")

#### SIMilarity Score ####
df_1 <- High_Res_IspH_AIIEIVDALDR_Spectra.RelativeAbundance
df_2 <- High_Res_IspH_AIIEIVDALDR_Spectra.RelativeAbundance
df_top10 <- AIEIVDQALDR_top10

SIM_function <- function(df_1, df_2, df_top10, name1, name2){
  
  
  df1_peptide <- peptide_ions_function(df_top10,df_1,name1) 
  
  df2_peptide <- peptide_ions_function(df_top10, df_2,name2)
  
  
  if(nrow(df1_peptide) >= nrow(df2_peptide)){
    spec1 = df1_peptide
    spec2 = df2_peptide
  }else{
    spec1 = df2_peptide
    spec2 = df1_peptide
  }
  
  sum_sqrt_Relative.Abundance = vector()
  #not_included = vector()
  
  for(i in 1:nrow(spec1)){
    
    print(i)
    
    if(is.na(trunc(spec2$m.z[i]))){
      print(paste0("i = ",i," Ion NOt here"))
      print(spec1$m.z[i])
      print(spec2$m.z[i])
    }
    
    if(length(which(trunc(spec1$m.z[i]) %in% trunc(spec2$m.z))) != 0){
      
      #spec1_intensity = spec1$Intensity[i]
      spec1_Relative.Abundance = spec1[,ncol(spec1)][i]
      
      #spec2_intensity = spec2$Intensity[which(trunc(spec2$m.z) == trunc(spec1$m.z[i]))]
      spec2_Relative.Abundance = spec2[,ncol(spec2)][which(trunc(spec2$m.z) == trunc(spec1$m.z[i]))]
      
      #combined_intensity = sqrt(spec1_intensity*spec2_intensity)
      sqrt_combined_Relative.Abundance = sqrt(spec1_Relative.Abundance*spec2_Relative.Abundance)
      
      #sqrt_intensity = append(sqrt_intensity, combined_intensity)
      sum_sqrt_Relative.Abundance = append(sum_sqrt_Relative.Abundance, sqrt_combined_Relative.Abundance)
      print(sum_sqrt_Relative.Abundance)
      
    } else{
      
      print(paste0("i = ",i," Ion NOt here"))
      
      spec1_Relative.Abundance = spec1[,ncol(spec1)][i]
      
      spec2_Relative.Abundance = 0
      
      combined_Relative.Abundance = sqrt(spec1_Relative.Abundance*spec2_Relative.Abundance)
      
      sum_sqrt_Relative.Abundance = append(sum_sqrt_Relative.Abundance, combined_Relative.Abundance)
      print(sum_sqrt_Relative.Abundance)
      
    }
  }
  
  I_spec1 <- intensity_function(df_top10, df1_peptide)
  I_spec2 <- intensity_function(df_top10, df2_peptide)
  
  numerator = sum(sum_sqrt_Relative.Abundance)
  denomenator = sqrt(sum(I_spec1)*sum(I_spec2))
  
  SIM = numerator/denomenator
  
  return(SIM)
  
}

SIM_function(High_Res_IspH_AIIEIVDALDR_Spectra.RelativeAbundance,High_Res_IspH_AIIEIVDALDR_Spectra.RelativeAbundance, AIEIVDQALDR_top10,"High_Res_IspH_AIIEIVDALDR_Spectra.RelativeAbundance","High_Res_IspH_AIIEIVDALDR_Spectra.RelativeAbundance")
SIM_function(High_Res_IspH_AIIEIVDALDR_Spectra.RelativeAbundance,Low_Res_IspH_AIIEIVDALDR_Spectra.RelaiveAbundance,AIEIVDQALDR_top10,"High_Res_IspH_AIIEIVDALDR_Spectra.RelativeAbundance","Low_Res_IspH_AIIEIVDALDR_Spectra.RelaiveAbundance")
SIM_function(High_Res_IspH_AIIEIVDALDR_Spectra.RelativeAbundance,FAIMS_AIIEIVDALDR_Spectra_RelativeAbundance, AIEIVDQALDR_top10,"High_Res_IspH_AIIEIVDALDR_Spectra.RelativeAbundance","FAIMS_AIIEIVDALDR_Spectra_RelativeAbundance")
SIM_function(High_Res_IspH_AIIEIVDALDR_Spectra.RelativeAbundance,noFAIMS_AIIEIVDALDR_Spectra_RelativeAbundance,AIEIVDQALDR_top10,"High_Res_IspH_AIIEIVDALDR_Spectra.RelativeAbundance","noFAIMS_AIIEIVDALDR_Spectra_RelativeAbundance")

df_matrix <- matrix(nrow = 5, ncol = 5)
df_matrix_ETDIGVTGGGGQGK <- matrix(nrow = 5, ncol = 5)

for(i in 1:length(df.list)){
  df_n <- df.list[[i]]
  df_name <- names(df.list[i])
  for (j in 1:length(df.list)) {
    SIM_score = SIM_function(df_n,df.list[[j]], AIEIVDQALDR_top10, df_name, names(df.list[j]))
    df_matrix[i,j] <- SIM_score
  }
}

for(i in 1:length(df.list_ETDIGVTGGGGQGK)){
  df_n <- df.list_ETDIGVTGGGGQGK[[i]]
  df_name <- names(df.list_ETDIGVTGGGGQGK[i])
  for (j in 1:length(df.list_ETDIGVTGGGGQGK)) {
    SIM_score = SIM_function(df_n,df.list_ETDIGVTGGGGQGK[[j]], ETDIGVTGGGQGK_top10,df_name, names(df.list_ETDIGVTGGGGQGK[j]))
    df_matrix_ETDIGVTGGGGQGK[i,j] <- SIM_score
  }
}

colnames(df_matrix) <- c("PROSIT","High Res", "Low Res", "Infusion FAIMS", "Infusion")
rownames(df_matrix) <- c("PROSIT","High Res", "Low Res", "Infusion FAIMS", "Infusion")

colnames(df_matrix_ETDIGVTGGGGQGK) <- c("PROSIT","High Res", "Low Res", "Infusion FAIMS", "Infusion")
rownames(df_matrix_ETDIGVTGGGGQGK) <- c("PROSIT","High Res", "Low Res", "Infusion FAIMS", "Infusion")

library(RColorBrewer)
#cols <- brewer.pal(4,"BrBG")
scaleRYG <- colorRampPalette(c("#F3A5BF","#15688E"), space = "rgb")(100)



pdf("F:/Projects/Proteomics/Zymomona/FAIMS/Figures/FromR/SimilarityScore_AIEIVDQALDR_ETDIGVTGGGQGK_v3.pdf",height = 10, width = 15)
par(mfrow = c(1,2))
corrplot(df_matrix, type = "lower", method = "color", #col = scaleRYG, 
         addCoef.col = "white",
         tl.col = "black", tl.srt = 45,
         cl.lim = c(0,1))

corrplot(df_matrix_ETDIGVTGGGGQGK, type = "lower", method = "color", #col = scaleRYG, 
         addCoef.col = "white",
         tl.col = "black", tl.srt = 45,
         cl.lim = c(0,1))
dev.off()


##### Spectral Similarity Function from Jesse ######

library(OrgMassSpecR)
 
df1_peptide <- peptide_ions_function(AIEIVDQALDR_top10,High_Res_IspH_AIIEIVDALDR_Spectra.RelativeAbundance,"High_Res_IspH_AIIEIVDALDR_Spectra.RelativeAbundance") 

df2_peptide <- peptide_ions_function(AIEIVDQALDR_top10, noFAIMS_AIIEIVDALDR_Spectra_RelativeAbundance,"noFAIMS_AIIEIVDALDR_Spectra_RelativeAbundance")

spectrumSimilarity_plotting_dp_function <- function(df_top10, df_1, df_2, name1, name2){
  
  df1_peptide <- peptide_ions_function(df_top10,df_1,name1)
  df2_peptide <- peptide_ions_function(df_top10, df_2,name2)
  
  dp <- SpectrumSimilarity(df1_peptide, df2_peptide,
                           top.label = name1,
                           bottom.label = name2)
  
  return <- dp
  
}

dotproduct = spectrumSimilarity_plotting_dp_function(AIEIVDQALDR_top10,PROSIT_IspH_AIIEIVDALDR_Spectra,High_Res_IspH_AIIEIVDALDR_Spectra.RelativeAbundance,
                                        "PROSIT_IspH_AIIEIVDALDR_Spectra","High_Res_IspH_AIIEIVDALDR_Spectra.RelativeAbundance")

j_matrix <- matrix(nrow = 5, ncol = 5)

pdf("F:/Projects/Proteomics/Zymomona/FAIMS/Figures/FromR/Jesse_Similarity_Plots.pdf",height = 10, width = 15)
par(mfrow = c(1,2))
for(i in 1:length(df.list)){
  df_n <- df.list[[i]]
  df_name <- names(df.list[i])
  for (j in 1:length(df.list)) {
    dotproduct_v2 = spectrumSimilarity_plotting_dp_function(AIEIVDQALDR_top10, df_n, df.list[[j]], df_name, names(df.list[j]))
    j_matrix[i,j] <- dotproduct_v2
  }
}
dev.off()

j_matrix_ETDIGVTGGGGQGK <- matrix(nrow = 5, ncol = 5)
pdf("F:/Projects/Proteomics/Zymomona/FAIMS/Figures/FromR/Jesse_Similarity_Plots_ETDIGVTGGGQGK.pdf",height = 10, width = 15)
par(mfrow = c(1,2))
for(i in 1:length(df.list)){
  df_n <- df.list_ETDIGVTGGGGQGK[[i]]
  df_name <- names(df.list_ETDIGVTGGGGQGK[i])
  for (j in 1:length(df.list_ETDIGVTGGGGQGK)) {
    dotproduct_v2 = spectrumSimilarity_plotting_dp_function(ETDIGVTGGGQGK_top10, df_n, df.list_ETDIGVTGGGGQGK[[j]], df_name, names(df.list[j]))
    j_matrix_ETDIGVTGGGGQGK[i,j] <- dotproduct_v2
  }
}
dev.off()

colnames(j_matrix) <- c("PROSIT","High Res", "Low Res", "Infusion FAIMS", "Infusion")
rownames(j_matrix) <- c("PROSIT","High Res", "Low Res", "Infusion FAIMS", "Infusion")

colnames(j_matrix_ETDIGVTGGGGQGK) <- c("PROSIT","High Res", "Low Res", "Infusion FAIMS", "Infusion")
rownames(j_matrix_ETDIGVTGGGGQGK) <- c("PROSIT","High Res", "Low Res", "Infusion FAIMS", "Infusion")


pdf("F:/Projects/Proteomics/Zymomona/FAIMS/Figures/FromR/Jesses_AIEIVDQALDR_ETDIGVTGGGQGK_heatmaps.pdf",height = 10, width = 15)
par(mfrow = c(1,2))
corrplot(j_matrix, type = "lower", method = "color", #col = scaleRYG, 
         addCoef.col = "white",
         tl.col = "black", tl.srt = 45)

corrplot(j_matrix_ETDIGVTGGGGQGK, type = "lower", method = "color", #col = scaleRYG, 
         addCoef.col = "white",
         tl.col = "black", tl.srt = 45)
dev.off()
