library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(corrplot)
library(ggcorrplot)


##### Library Comparison ####
ProteinProspector_AIEIVQALDR <- read.csv("H:/Projects/Proteomics/Zymomona/FAIMS/DataAnalysis/ProteinProspector_AIEIVDQALDR.csv", 
                                         header = TRUE, sep = ",", stringsAsFactors = FALSE)
ProteinProspector_ETDIGVTGGGQGK <- read.csv("H:/Projects/Proteomics/Zymomona/FAIMS/DataAnalysis/ProteinProspector_ETDIGVTGGGQGK.csv", 
                                            header = TRUE, sep = ",", stringsAsFactors = FALSE)

#Library Spectra
#PROSIT
PROSIT_IspH_AIIEIVDALDR_Spectra <- read.csv("H:/Projects/Proteomics/Zymomona/FAIMS/DataAnalysis/LibraryComparison/PROSIT_LibrarySpectra_AIEIVDQALDR.csv", 
                                            header = TRUE, sep = ",", stringsAsFactors = FALSE)
colnames(PROSIT_IspH_AIIEIVDALDR_Spectra) <- c("m.z", "Relative.Abundance")

PROSIT_IspG_ETDIGVTGGGGQGK_Spectra <- read.csv("H:/Projects/Proteomics/Zymomona/FAIMS/DataAnalysis/LibraryComparison/PROSIT_LibrarySpectra_ETDIGVTGGGGQGK.csv", 
                                               header = TRUE, sep = ",", stringsAsFactors = FALSE)
colnames(PROSIT_IspG_ETDIGVTGGGGQGK_Spectra) <- c("m.z", "Relative.Abundance")

#### HIGH_Res data upload ####
High_Res_IspH_AIIEIVDALDR_Spectra <- read.csv("H:/Projects/Proteomics/Zymomona/FAIMS/DataAnalysis/LibraryComparison/HighRes_Intensity_LibrarySpectra_AIEIVDQALDR.csv", 
                                              header = TRUE, sep = ",", stringsAsFactors = FALSE)
High_Res_IspH_AIIEIVDALDR_Spectra.RelativeAbundance <- High_Res_IspH_AIIEIVDALDR_Spectra
High_Res_IspH_AIIEIVDALDR_Spectra.RelativeAbundance$Relative.Abundance <- High_Res_IspH_AIIEIVDALDR_Spectra.RelativeAbundance$Intensity/(max(High_Res_IspH_AIIEIVDALDR_Spectra.RelativeAbundance$Intensity))
colnames(High_Res_IspH_AIIEIVDALDR_Spectra.RelativeAbundance) <- c("m.z", "Intensity", "Relative.Abundance")



High_Res_IspG_ETDIGVTGGGGQGK_Spectra <- read.csv("H:/Projects/Proteomics/Zymomona/FAIMS/DataAnalysis/LibraryComparison/HighRes_LibrarySpectra_ETDIGVTGGGGQGK.csv", 
                                                 header = TRUE, sep = ",", stringsAsFactors = FALSE)
High_Res_IspG_ETDIGVTGGGGQGK_Spectra.RelativeAbundance = High_Res_IspG_ETDIGVTGGGGQGK_Spectra
High_Res_IspG_ETDIGVTGGGGQGK_Spectra.RelativeAbundance$Relative.Abundance <- High_Res_IspG_ETDIGVTGGGGQGK_Spectra.RelativeAbundance$Intensity/(max(High_Res_IspG_ETDIGVTGGGGQGK_Spectra.RelativeAbundance$Intensity))
colnames(High_Res_IspG_ETDIGVTGGGGQGK_Spectra.RelativeAbundance) <- c("m.z", "Intensity", "Relative.Abundance")



#### LOW_Res data upload ####
Low_Res_IspH_AIIEIVDALDR_Spectra <- read.csv("H:/Projects/Proteomics/Zymomona/FAIMS/DataAnalysis/LibraryComparison/LowRes_LibrarySpectra_AIEIVDQALDR.csv", 
                                             header = TRUE, sep = ",", stringsAsFactors = FALSE)
Low_Res_IspH_AIIEIVDALDR_Spectra.RelaiveAbundance <- Low_Res_IspH_AIIEIVDALDR_Spectra
Low_Res_IspH_AIIEIVDALDR_Spectra.RelaiveAbundance$Relative.Abundance <- Low_Res_IspH_AIIEIVDALDR_Spectra.RelaiveAbundance$Intensity/(max(Low_Res_IspH_AIIEIVDALDR_Spectra.RelaiveAbundance$Intensity))
colnames(Low_Res_IspH_AIIEIVDALDR_Spectra.RelaiveAbundance) <- c("m.z", "Intensity", "Relative.Abundance")


Low_Res_IspG_ETDIGVTGGGGQGK_Spectra <- read.csv("H:/Projects/Proteomics/Zymomona/FAIMS/DataAnalysis/LibraryComparison/LowRes_LibrarySpectra_ETDIGVTGGGGQGK.csv", 
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


##################################################



# df_library = ProteinProspector_ETDIGVTGGGQGK
# df = FAIMS_ETDIGVTGGGQGK_Spectra_RelativeAbundance
spectra_plotting <- function(df_library, df){
  # extract the library peptide ions into object fragment
  fragment <- numeric()
  for (i in 1:nrow(df_library)){
    if(df_library[i,2] != 0 ){
      fragment <- append(fragment,as.numeric(df_library[i,1]))
    }
  }
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
    }
  
  # extract the m/z values using the index were pattern matching was found
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

##### Plot spectra #####
prosit_AIEIVDQALDR <- spectra_plotting(ProteinProspector_AIEIVQALDR, PROSIT_IspH_AIIEIVDALDR_Spectra)
PA <- prosit_AIEIVDQALDR + labs(title = "Prosit Theoretical Spectra - AIEVDQALDR", y = "Relative Abundance", x = "m/z")

high_res_AIEIDVDQDALDR <- spectra_plotting(ProteinProspector_AIEIVQALDR, High_Res_IspH_AIIEIVDALDR_Spectra.RelativeAbundance)
HA <- high_res_AIEIDVDQDALDR + labs(title = "High Resolution AIEIVDQDALDR Spectra", y = "Relative Abundance", x ="m/z")

low_res_AIEIDVDQDALDR <- spectra_plotting(ProteinProspector_AIEIVQALDR, Low_Res_IspH_AIIEIVDALDR_Spectra.RelaiveAbundance)
LA <- low_res_AIEIDVDQDALDR + labs(title = "Low Resolution AIEIVDQDALDR Spectra", y = "Relative Abundance", x ="m/z")
#labs(title = "Low Resolution ETDIGVTGGGGQGK Spectra", y = "Relative Abundance", x ="m/z")

prosit_ETDIGVTGGGGQGK <- spectra_plotting(ProteinProspector_ETDIGVTGGGQGK, PROSIT_IspG_ETDIGVTGGGGQGK_Spectra)
PE <- prosit_ETDIGVTGGGGQGK + labs(title = "Prosit Theoretical Spectra - ETDIGVTGGGQGK", y= "Relative Abundance", x="m/z")

high_res_ETDIGVTGGGGQGK <- spectra_plotting(ProteinProspector_ETDIGVTGGGQGK, High_Res_IspG_ETDIGVTGGGGQGK_Spectra.RelativeAbundance)
HE <- high_res_ETDIGVTGGGGQGK + labs(title = "High Resolution ETDIGVTGGGGQGK Spectra", y = "Relative Abundance", x ="m/z")

low_res_ETDIGVTGGGGQGK <- spectra_plotting(PROSIT_IspG_ETDIGVTGGGGQGK_Spectra, Low_Res_IspG_ETDIGVTGGGGQGK_Spectra.RelativeAbundance)
LE <- low_res_ETDIGVTGGGGQGK + labs(title = "Low Resolution ETDIGVTGGGGQGK Spectra", y = "Relative Abundance", x ="m/z")


##### Infusion Comparison to Library ######
# The following functin infusion_spectra_prep is meant to format the infusion spectra for generating 
# the spectral plots

# Infusion_df <- noFAIMS_AIIEIVDALDR_Spectra
# max_Intensity <- 200000

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
#FAIMS
FAIMS_AIIEIVDALDR_Spectra <- read.csv("H:/Projects/Proteomics/Zymomona/FAIMS/DataAnalysis/LibraryComparison/20200126_ZymoFAIMS_IW2_R500K_MI502_AGC1e06_2.csv", 
                                      header = TRUE, sep = ",", stringsAsFactors = FALSE)

FAIMS_AIIEIVDALDR_Spectra_RelativeAbundance <- infusion_spectra_prep(FAIMS_AIIEIVDALDR_Spectra, 200000)

FAIMS_AIEIDVDQDALDR <- spectra_plotting(ProteinProspector_AIEIVQALDR, FAIMS_AIIEIVDALDR_Spectra_RelativeAbundance)
A1 <- FAIMS_AIEIDVDQDALDR + ylim(0,1) + labs(title = "FAIMS AIEIVDQDALDR Spectra", y = "Relative Abundance", x ="m/z")

#with out FAIMS
noFAIMS_AIIEIVDALDR_Spectra <- read.csv("H:/Projects/Proteomics/Zymomona/FAIMS/DataAnalysis/LibraryComparison/20200125_Zymo_IW2_R500K_MI502_AGC1e06_1.csv", 
                                        header = TRUE, sep = ",", stringsAsFactors = FALSE)

noFAIMS_AIIEIVDALDR_Spectra_RelativeAbundance <- infusion_spectra_prep(noFAIMS_AIIEIVDALDR_Spectra, 200000)

noFAIMS_AIEIDVDQDALDR <- spectra_plotting(ProteinProspector_AIEIVQALDR, noFAIMS_AIIEIVDALDR_Spectra_RelativeAbundance)
A2 <- noFAIMS_AIEIDVDQDALDR + ylim(0,1) + labs(title = "no FAIMS AIEIVDQDALDR Spectra", y = "Relative Abundance", x ="m/z")

# ETDIGVTGGGQGK #

#FAIMS
FAIMS_ETDIGVTGGGQGK_Spectra <- read.csv("H:/Projects/Proteomics/Zymomona/FAIMS/DataAnalysis/LibraryComparison/20200126_ZymoFAIMS_IW2_R500K_MI502_AGC1e06_2_ETDIGVTGGGQGK.csv", 
                                        header = TRUE, sep = ",", stringsAsFactors = FALSE)

FAIMS_ETDIGVTGGGQGK_Spectra_RelativeAbundance <- infusion_spectra_prep(FAIMS_ETDIGVTGGGQGK_Spectra, 200000)

FAIMS_ETDIGVTGGGQGK <- spectra_plotting(ProteinProspector_ETDIGVTGGGQGK, FAIMS_ETDIGVTGGGQGK_Spectra_RelativeAbundance)
E1 <- FAIMS_ETDIGVTGGGQGK + ylim(0,1) + labs(title = "FAIMS ETDIGVTGGGQGK Spectra", y = "Relative Abundance", x ="m/z")

#with out FAIMS
noFAIMS_ETDIGVTGGGQGK_Spectra <- read.csv("H:/Projects/Proteomics/Zymomona/FAIMS/DataAnalysis/LibraryComparison/20200125_Zymo_IW2_R500K_MI502_AGC1e06_1_ETDIGVTGGGQGK.csv", 
                                          header = TRUE, sep = ",", stringsAsFactors = FALSE)

noFAIMS_ETDIGVTGGGQGK_Spectra_RelativeAbundance <- infusion_spectra_prep(noFAIMS_ETDIGVTGGGQGK_Spectra, 200000)

noFAIMS_ETDIGVTGGGGQGK <- spectra_plotting(PROSIT_IspG_ETDIGVTGGGGQGK_Spectra, noFAIMS_ETDIGVTGGGQGK_Spectra_RelativeAbundance)
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


pdf("H:/Projects/Proteomics/Zymomona/FAIMS/Figures/FromR/SpectraComparisons_PROSIT_HIGHRES_LOWRES.pdf")
grid.arrange(PAgrob, PEgrob, HAgrob, HEgrob, LAgrob, LEgrob, ncol =2)
dev.off()


ggsave("H:/Projects/Proteomics/Zymomona/FAIMS/Figures/FromR/SpectraComparisons_PROSIT_HIGHRES_LOWRES.pdf",
       grid.arrange(PAgrob, PEgrob, HAgrob, HEgrob, LAgrob, LEgrob, ncol =2),
       width = 10,
       height = 11.5,
       units = "in")

ggsave("H:/Projects/Proteomics/Zymomona/FAIMS/Figures/FromR/SpectraComparisons_FAIMS_noFAIMS.pdf",
       grid.arrange(A1grob, E1grob, A2grob, E2grob, ncol =2),
       width = 10,
       height = 8.5,
       units = "in")

ggsave("H:/Projects/Proteomics/Zymomona/FAIMS/Figures/FromR/SpectraComparisons_all.pdf",
       grid.arrange(PAgrob, PEgrob, HAgrob, HEgrob, LAgrob, LEgrob, A1grob, E1grob, A2grob, E2grob, ncol =2),
       width = 10,
       height = 14.5,
       units = "in")

#### TIC Explained ######
# To explain the TIC associated to the peptide fragments the function tic_explained_function

# df_proteinProspect <- ProteinProspector_AIEIVQALDR
# df <- FAIMS_AIIEIVDALDR_Spectra_RelativeAbundance

tic_explained_function <- function(df_proteinProspect, df){
  
  transitions <- numeric()
  
  df.mass.to.charge.char <- as.character(trunc(df_proteinProspect$mass.to.charge))
  
  
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

Tic_explained_High_Res_AIEIDVQALDR <- tic_explained_function(ProteinProspector_AIEIVQALDR, High_Res_IspH_AIIEIVDALDR_Spectra)

Tic_explained_Low_Res_AIEIDVQALDR <- tic_explained_function(ProteinProspector_AIEIVQALDR,Low_Res_IspH_AIIEIVDALDR_Spectra)

Tic_explained_High_Res_ETDIGVTGGGQGK <- tic_explained_function(ProteinProspector_ETDIGVTGGGQGK, High_Res_IspG_ETDIGVTGGGGQGK_Spectra)

Tic_explained_Low_Res_ETDIGVTGGGQGK <- tic_explained_function(ProteinProspector_ETDIGVTGGGQGK, Low_Res_IspG_ETDIGVTGGGGQGK_Spectra)

### AIEIVDQALDR ###
colnames(FAIMS_AIIEIVDALDR_Spectra_RelativeAbundance) <- c("m.z", "Intensity", "Relative.Abundance")
Tic_explained_FAIMS_AIEIDVQALDR <- tic_explained_function(ProteinProspector_AIEIVQALDR,FAIMS_AIIEIVDALDR_Spectra_RelativeAbundance)

colnames(noFAIMS_AIIEIVDALDR_Spectra_RelativeAbundance) <- c("m.z", "Intensity", "Relative.Abundance")
Tic_explained_noFAIMS_AIEIDVQALDR <- tic_explained_function(ProteinProspector_AIEIVQALDR, noFAIMS_AIIEIVDALDR_Spectra_RelativeAbundance)

### ETDIGVTGGGQGK ###
colnames(FAIMS_ETDIGVTGGGQGK_Spectra_RelativeAbundance) <- c("m.z", "Intensity", "Relative.Abundance")
Tic_explained_FAIMS_ETDIGVTGGGQGK <- tic_explained_function(ProteinProspector_ETDIGVTGGGQGK,FAIMS_ETDIGVTGGGQGK_Spectra_RelativeAbundance)

colnames(noFAIMS_ETDIGVTGGGQGK_Spectra_RelativeAbundance) <- c("m.z", "Intensity", "Relative.Abundance")
Tic_explained_noFAIMS_ETDIGVTGGGQGK <- tic_explained_function(ProteinProspector_ETDIGVTGGGQGK, noFAIMS_ETDIGVTGGGQGK_Spectra_RelativeAbundance)



##### Dot Product ####
# df_proteinProspect <- ProteinProspector_AIEIVQALDR
# df <- High_Res_IspH_AIIEIVDALDR_Spectra

peptide_ions_function <- function(df_proteinProspect, df){
  
  transitions <- numeric()
  
  df.mass.to.charge.char <- as.character(trunc(df_proteinProspect$mass.to.charge))
  
  
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
      
      keep.intensity <- max(df$Intensity[which(df$m.z %in% mass.to.charge[comparison_intensity])])
      
      keep <- df$m.z[which(df$Intensity == keep.intensity)]
      
      final_list <- append(final_list, keep)
      #print(final_list)
    }
    
  }
  
  peptide <- df[which(df$m.z %in% final_list),]
  
  return(peptide)
  
}



# df_proteinProspect <- ProteinProspector_AIEIVQALDR
# df <- High_Res_IspH_AIIEIVDALDR_Spectra

intensity_function <- function(df_proteinProspect, df){
  
  transitions <- numeric()
  
  df.mass.to.charge.char <- as.character(trunc(df_proteinProspect$mass.to.charge))
  
  
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
      
      keep.intensity <- max(df$Intensity[which(df$m.z %in% mass.to.charge[comparison_intensity])])
      
      keep <- df$m.z[which(df$Intensity == keep.intensity)]
      
      final_list <- append(final_list, keep)
      #print(final_list)
    }
    
  }
  
  intensity_peptide <- df$Intensity[which(df$m.z %in% final_list)]
  #intensity_peptide_sum <- sum(intensity_peptide)
  
  return(intensity_peptide)
  
}

# df_1 <- High_Res_IspG_ETDIGVTGGGGQGK_Spectra
# df_2 <- noFAIMS_ETDIGVTGGGQGK_Spectra_RelativeAbundance
# df_ProteinProspector <- ProteinProspector_ETDIGVTGGGQGK

DP_function <- function(df_1, df_2, df_ProteinProspector){
  
  df1_peptide <- peptide_ions_function(df_ProteinProspector,df_1)
  df2_peptide <- peptide_ions_function(df_ProteinProspector, df_2)
  
  if(nrow(df1_peptide) >= nrow(df2_peptide)){
    spec1 = df1_peptide
    spec2 = df2_peptide
  }else{
    spec1 = df2_peptide
    spec2 = df1_peptide
  }
  
  sum_intensity_combined = vector()
  not_included = vector()
  
  for(i in 1:nrow(spec1)){
    
    print(i)
    
    if(is.na(trunc(spec2$m.z[i]))){
      print(paste0("i = ",i," Ion NOt here"))
      print(spec1$m.z[i])
      print(spec2$m.z[i])
    }
    
    if(length(which(trunc(spec1$m.z[i]) %in% trunc(spec2$m.z))) != 0){
      
      spec1_intensity = spec1$Intensity[i]
      
      spec2_intensity = spec2$Intensity[which(trunc(spec2$m.z) == trunc(spec1$m.z[i]))]
      
      combined_intensity = spec1_intensity*spec2_intensity
      
      sum_intensity_combined = append(sum_intensity_combined, combined_intensity)
      print(sum_intensity_combined)
      
    } else{
      
      print(paste0("i = ",i," Ion NOt here"))
      
    }
  }
  
  
  I_spec1 <- intensity_function(df_ProteinProspector, df_1)
  I_spec2 <- intensity_function(df_ProteinProspector, df_2)
  
  numerator <- sum(sum_intensity_combined)
  denomenator <- sum(I_spec1^2)*sum(I_spec2^2)
  
  DP <- numerator/sqrt(denomenator)
  
  return(DP)
  
}


DP_function(High_Res_IspH_AIIEIVDALDR_Spectra,noFAIMS_AIIEIVDALDR_Spectra_RelativeAbundance, ProteinProspector_AIEIVQALDR)
DP_function(High_Res_IspG_ETDIGVTGGGGQGK_Spectra,noFAIMS_ETDIGVTGGGQGK_Spectra_RelativeAbundance, ProteinProspector_ETDIGVTGGGQGK)


df.list <- list(High_Res_IspH_AIIEIVDALDR_Spectra, Low_Res_IspH_AIIEIVDALDR_Spectra, 
                FAIMS_AIIEIVDALDR_Spectra_RelativeAbundance, noFAIMS_AIIEIVDALDR_Spectra_RelativeAbundance)

df.list_ETDIGVTGGGGQGK <- list(High_Res_IspH_ETDIGVTGGGGQGK_Spectra, Low_Res_IspH_ETDIGVTGGGGQGK_Spectra, 
                FAIMS_ETDIGVTGGGQGK_Spectra_RelativeAbundance, noFAIMS_ETDIGVTGGGQGK_Spectra_RelativeAbundance)


colnames(df.list_ETDIGVTGGGGQGK[[1]]) <- c("m.z","Intensity")
colnames(df.list_ETDIGVTGGGGQGK[[2]]) <- c("m.z","Intensity")

df_matrix <- matrix(nrow = 4, ncol = 4)
df_matrix_ETDIGVTGGGGQGK <- matrix(nrow = 4, ncol = 4)

for(i in 1:length(df.list)){
  df_n <- df.list[[i]]
  for (j in 1:length(df.list)) {
    dot_product = DP_function(df_n,df.list[[j]],ProteinProspector_AIEIVQALDR)
    df_matrix[i,j] <- dot_product
  }
}

for(i in 1:length(df.list_ETDIGVTGGGGQGK)){
  df_n <- df.list_ETDIGVTGGGGQGK[[i]]
  for (j in 1:length(df.list_ETDIGVTGGGGQGK)) {
    dot_product = DP_function(df_n,df.list_ETDIGVTGGGGQGK[[j]], ProteinProspector_ETDIGVTGGGQGK)
    df_matrix_ETDIGVTGGGGQGK[i,j] <- dot_product
  }
}

colnames(df_matrix) <- c("High Res", "Low Res", "Infusion FAIMS", "Infusion")
rownames(df_matrix) <- c("High Res", "Low Res", "Infusion FAIMS", "Infusion")

colnames(df_matrix_ETDIGVTGGGGQGK) <- c("High Res", "Low Res", "Infusion FAIMS", "Infusion")
rownames(df_matrix_ETDIGVTGGGGQGK) <- c("High Res", "Low Res", "Infusion FAIMS", "Infusion")

library(RColorBrewer)
#cols <- brewer.pal(4,"BrBG")
scaleRYG <- colorRampPalette(c("#F3A5BF","#15688E"), space = "rgb")(100)



pdf("H:/Projects/Proteomics/Zymomona/FAIMS/Figures/FromR/dotProduct_AIEIVDQALDR_ETDIGVTGGGQGK_v2.pdf",height = 10, width = 15)
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
# df_1 <- High_Res_IspH_AIIEIVDALDR_Spectra
# df_2 <- High_Res_IspH_AIIEIVDALDR_Spectra

SIM_function <- function(df_1, df_2, df_ProteinProspector){
  
   df1_peptide <- peptide_ions_function(df_ProteinProspector,df_1)
   df2_peptide <- peptide_ions_function(df_ProteinProspector, df_2)
   
   if(nrow(df1_peptide) >= nrow(df2_peptide)){
     spec1 = df1_peptide
     spec2 = df2_peptide
   }else{
     spec1 = df2_peptide
     spec2 = df1_peptide
   }
   
  sqrt_intensity = vector()
  not_included = vector()
  
   for(i in 1:nrow(spec1)){
     
     print(i)
     
     if(is.na(trunc(spec2$m.z[i]))){
       print(paste0("i = ",i," Ion NOt here"))
       print(spec1$m.z[i])
       print(spec2$m.z[i])
     }
     
     if(length(which(trunc(spec1$m.z[i]) %in% trunc(spec2$m.z))) != 0){
       
      spec1_intensity = spec1$Intensity[i]
      
      spec2_intensity = spec2$Intensity[which(trunc(spec2$m.z) == trunc(spec1$m.z[i]))]
      
      combined_intensity = sqrt(spec1_intensity*spec2_intensity)
         
      sqrt_intensity = append(sqrt_intensity, combined_intensity)

     } else{
       
       print(paste0("i = ",i," Ion NOt here"))
       
     }
   }
  
  I_spec1 <- intensity_function(df_ProteinProspector, df_1)
  I_spec2 <- intensity_function(df_ProteinProspector, df_2)
  
  numerator = sum(sqrt_intensity)
  denomenator = sqrt(sum(I_spec1)*sum(I_spec2))
  
  SIM = numerator/denomenator
  
  return(SIM)
  
}

SIM_function(High_Res_IspH_AIIEIVDALDR_Spectra,High_Res_IspH_AIIEIVDALDR_Spectra, ProteinProspector_AIEIVQALDR)
SIM_function(High_Res_IspH_AIIEIVDALDR_Spectra,Low_Res_IspH_AIIEIVDALDR_Spectra,ProteinProspector_AIEIVQALDR)
SIM_function(High_Res_IspH_AIIEIVDALDR_Spectra,FAIMS_AIIEIVDALDR_Spectra_RelativeAbundance, ProteinProspector_AIEIVQALDR)
SIM_function(High_Res_IspH_AIIEIVDALDR_Spectra,noFAIMS_AIIEIVDALDR_Spectra_RelativeAbundance)

df_matrix <- matrix(nrow = 4, ncol = 4)
df_matrix_ETDIGVTGGGGQGK <- matrix(nrow = 4, ncol = 4)

for(i in 1:length(df.list)){
  df_n <- df.list[[i]]
  for (j in 1:length(df.list)) {
    SIM_score = SIM_function(df_n,df.list[[j]], ProteinProspector_AIEIVQALDR)
    df_matrix[i,j] <- SIM_score
  }
}

for(i in 1:length(df.list_ETDIGVTGGGGQGK)){
  df_n <- df.list_ETDIGVTGGGGQGK[[i]]
  for (j in 1:length(df.list_ETDIGVTGGGGQGK)) {
    SIM_score = SIM_function(df_n,df.list_ETDIGVTGGGGQGK[[j]], ProteinProspector_ETDIGVTGGGQGK)
    df_matrix_ETDIGVTGGGGQGK[i,j] <- SIM_score
  }
}

colnames(df_matrix) <- c("High Res", "Low Res", "Infusion FAIMS", "Infusion")
rownames(df_matrix) <- c("High Res", "Low Res", "Infusion FAIMS", "Infusion")

colnames(df_matrix_ETDIGVTGGGGQGK) <- c("High Res", "Low Res", "Infusion FAIMS", "Infusion")
rownames(df_matrix_ETDIGVTGGGGQGK) <- c("High Res", "Low Res", "Infusion FAIMS", "Infusion")

library(RColorBrewer)
#cols <- brewer.pal(4,"BrBG")
scaleRYG <- colorRampPalette(c("#F3A5BF","#15688E"), space = "rgb")(100)



pdf("H:/Projects/Proteomics/Zymomona/FAIMS/Figures/FromR/SimilarityScore_AIEIVDQALDR_ETDIGVTGGGQGK_v2.pdf",height = 10, width = 15)
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
