library(ggplot2)
library(plotly)
library(pheatmap)
library(grid)
library(gridExtra)
library(lattice)

#### WT_3c_FAIMS_ISPHpeptides_infusion_20191112 ####

WT_FAIMS_ISPH <- read.csv("I:/EAT_BackUP/Trujillo_34362342_44636349_SN_Data/P2/Complete/Projects/Proteomics/Zymomona/WT_3c_FAIMS_ISPHpeptides_infusion_20191112.csv", 
                                            header = TRUE, sep = ",", stringsAsFactors = FALSE)

WT_FAIMS_ISPH$Peptide <- as.factor(WT_FAIMS_ISPH$Peptide)

p <- ggplot(WT_FAIMS_ISPH, aes(x=Retention.Time, y=Intensity, group = Peptide)) +
  geom_line(aes(color=Peptide), size=2) + 
  scale_color_manual(values=c('#BB67BC', '#8D7ADD', '#68D6B9', '#54ADD3'))  
p




#### ISPH_3a_FAIMS_ISPHpeptides_infusion_20191112 ####

# ISPH_FAIMS_ISPH <- read.csv("I:/EAT_BackUP/Trujillo_34362342_44636349_SN_Data/P2/Complete/Projects/Proteomics/Zymomona/ISPH_3a_FAIMS_ISPHpeptides_infusion_20191112.csv", 
#                           header = TRUE, sep = ",", stringsAsFactors = FALSE)

ISPH_FAIMS_ISPH <- read.csv("I:/EAT_BackUP/Trujillo_34362342_44636349_SN_Data/P2/Complete/Projects/Proteomics/Zymomona/ISPH_3a_FAIMS_infusioon_240K_Is_04iso_peptide.csv",
                            header = TRUE, sep = ",", stringsAsFactors = FALSE)

ISPH_plot <- function(x){
  x[,1] <- as.factor(x[,1])
  p <- ggplot(x, aes(x=Retention.Time, y=Intensity, group = Peptide)) +
    geom_line(aes(color=Peptide), size=2) + 
    scale_y_continuous(breaks = seq(0,60000,10000), limit = c(0,60000), labels = function(x)  format(x, scientific = TRUE)) + #max(x[,3])) +
    xlim(0.9, 2.1)+
    scale_color_manual(values=c('#BB67BC', '#8D7ADD', '#68D6B9', '#54ADD3'))+
    facet_wrap(~Replicate,nrow = 1)+
    theme_classic()
  return(p)
}

ISPG_plot <- function(x){
  x[,1] <- as.factor(x[,1])
  p <- ggplot(x, aes(x=Retention.Time, y=Intensity, group = Peptide)) +
    geom_line(aes(color=Peptide), size=2) + 
    scale_y_continuous(breaks = seq(0,600000,100000), limit = c(0,600000)) + #max(x[,3])) +
    xlim(0.9, 2.1)+
    scale_color_manual(values=c('#BB67BC', '#8D7ADD', '#68D6B9', '#54ADD3'))+
    facet_wrap(~Replicate,nrow = 1)+
    theme_classic()
  return(p)
}

ISPH_plot(ISPH_FAIMS_ISPH)

x <- ISPH_FAIMS_ISPH
x[,1] <- as.factor(x[,1])
ggplot(x, aes(x=Retention.Time, y=Intensity, group = Peptide)) +
  geom_line(aes(color=Peptide), size=2) + 
  #scale_y_continuous(breaks = seq(0,60000,10000), limit = c(0,60000), labels = function(x)  format(x, scientific = TRUE)) + #max(x[,3])) +
  xlim(1.0, 2.0)+
  scale_color_manual(values=c('#BB67BC', '#8D7ADD', '#68D6B9', '#54ADD3')) +
  theme_classic()+
  labs(title = "1 min tryptic digest infusion with FAIMS-PRM" ,subtitle = "Strain over-expressinig IspH", y = "Intensity", x ="Time (min)")

##### 8 peptide infusion ####

EightPeptide_3rep_Infusion <- read.csv("F:/Projects/Proteomics/Zymomona/FAIMS/DataAnalysis/FAIMS_8peptides_SpikedMatrix_Rep123_62_5fmol_edit.csv",
                                  header = TRUE, sep = ",", stringsAsFactors = FALSE)
EightPeptide_3rep_Infusion$Peptide <- as.factor(EightPeptide_Infusion$Peptide)
EightPeptide_3rep_Infusion$Replicate <- as.factor(EightPeptide_3rep_Infusion$Replicate)



ggplot(EightPeptide_3rep_Infusion, aes(x=Retention.Time, y=Intensity, group = Peptide)) +
  geom_point(aes(shape=Peptide, color=Replicate), size = 2)+
  #geom_line(aes(color=Peptide), size=1) + 
  scale_y_continuous(breaks = seq(0,3000000,500000), limit = c(0,3000000), labels = function(x)  format(x, scientific = TRUE)) +
  xlim(1.0, 2.0)+
  theme_classic()+
  labs(title = "8 peptides",x = "Time (min)")



EightPeptide_Infusion <- read.csv("F:/Projects/Proteomics/Zymomona/FAIMS/DataAnalysis/FAIMS_8peptides_SpikedMatrix_Rep3_62_5fmol_edit.csv",
                            header = TRUE, sep = ",", stringsAsFactors = FALSE)
 EightPeptide_Infusion$Peptide <- as.factor(EightPeptide_Infusion$Peptide)

p_8infusions <- ggplot(EightPeptide_Infusion, aes(x=Retention.Time, y=Intensity, group = Peptide)) +
  geom_point(aes(shape=Peptide, color=Peptide), size = 2)+
  geom_line(aes(color=Peptide), size=1) + 
  scale_y_continuous(breaks = seq(0,3000000,500000), limit = c(0,3000000), labels = function(x)  format(x, scientific = TRUE)) +
  xlim(1.0, 2.0)+
  theme_classic()+
  labs(title = "8 peptides",x = "Time (min)")
p_8infusions

FAIMS_8peptides_ETD_10transitions <- read.csv("F:/Projects/Proteomics/Zymomona/FAIMS/DataAnalysis/FAIMS_8peptides_SpikedMatrix_Rep3_62_5fmol_Transitions_ETDIGVTGGGQGK.csv", 
                                                 header = TRUE, sep = ",", stringsAsFactors = FALSE)
colnames(FAIMS_8peptides_ETD_10transitions) <- c("Time","Intensity","Transition")
x <- FAIMS_8peptides_ETD_10transitions
x[,3] <- as.factor(x[,3])
p_ETD_transitions  <- ggplot(x, aes(x=Time, y=Intensity, group = Transition)) +
  geom_line(aes(color=Transition), size=1) + 
  scale_y_continuous(breaks = seq(0,7050000,150000), limit = c(0,750000), labels = function(x)  format(x, scientific = TRUE)) +
  xlim(1.0, 2.0)+
  theme_classic()+
  labs(title = "ETDIG peptide transitions" , y = "Intensity", x ="Time (min)")

p_8infusions_grob <- ggplotGrob(p_8infusions)
p_ETD_transitions_grob <- ggplotGrob(p_ETD_transitions)

pdf("F:/Projects/Proteomics/Zymomona/FAIMS/Figures/FromR/SpectraComparisons_PROSIT_HIGHRES_LOWRES_IW6.pdf")
grid.draw(rbind(p_8infusions_grob, p_ETD_transitions_grob, size = "last"))
dev.off()

#### ISPH_3a_FAIMS_ISPHpeptides_infusion_240K ####

ISPH_FAIMS_240K <- read.csv("I:/EAT_BackUP/Trujillo_34362342_44636349_SN_Data/P2/Complete/Projects/Proteomics/Zymomona/ISPH_3a_FAIMS_ISPHpeptides_infusion_240K.csv", 
                            header = TRUE, sep = ",", stringsAsFactors = FALSE)

ISPH_plot(ISPH_FAIMS_240K)

##### Analysis Time #####

Analysis_Time <- read.csv("F:/Projects/Proteomics/Zymomona/FAIMS/DataAnalysis/AnalysisTime_LC_INFUSION.csv", 
                           header = TRUE, sep = ",", stringsAsFactors = FALSE)

ggplot(Analysis_Time, aes(x=Analysis, y=Time)) +
  geom_bar(stat="identity", width = 0.5)+
  theme_classic()
  

##### Infusion Spectra extraction of AIEVDQALDR #####

# FAIMS_IspH_AIIEIVDALDR_10transitions <- read.csv("I:/EAT_BackUP/Trujillo_34362342_44636349_SN_Data/P2/Complete/Projects/Proteomics/Zymomona/ISPH_3a_FAIMS_ISPHpeptides_infusion_LongTransitions_20191112.csv", 
#                                                 header = TRUE, sep = ",", stringsAsFactors = FALSE)

FAIMS_IspH_AIIEIVDALDR_10transitions <- read.csv("I:/EAT_BackUP/Trujillo_34362342_44636349_SN_Data/P2/Complete/Projects/Proteomics/Zymomona/ISPH_3a_FAIMS_infusioon_240K_Is_04iso_LongTransition.csv", 
                                                 header = TRUE, sep = ",", stringsAsFactors = FALSE)
colnames(FAIMS_IspH_AIIEIVDALDR_10transitions) <- c("Time","Intensity","Transition")
x <- FAIMS_IspH_AIIEIVDALDR_10transitions
x[,3] <- as.factor(x[,3])
ggplot(x, aes(x=Time, y=Intensity, group = Transition)) +
  geom_line(aes(color=Transition), size=2) + 
  #scale_y_continuous(breaks = seq(0,60000,10000), limit = c(0,60000), labels = function(x)  format(x, scientific = TRUE)) + #max(x[,3])) +
  xlim(1.0, 2.0)+
  #scale_color_manual(values=c('#BB67BC', '#8D7ADD', '#68D6B9', '#54ADD3')) +
  theme_classic()+
  labs(title = "1 min tryptic digest infusion with FAIMS-PRM" ,subtitle = "Strain over-expressinig IspH", y = "Intensity", x ="Time (min)")


# FAIMS_IspH_AIIEIVDALDR_Spectra <- read.csv("I:/EAT_BackUP/Trujillo_34362342_44636349_SN_Data/P2/Complete/Projects/Proteomics/Zymomona/ISPH_3a_FAIMS_ISPHpeptides_infusion_spectra_allTransitions_20191112.csv", 
#                                                  header = TRUE, sep = ",", stringsAsFactors = FALSE)

FAIMS_IspH_AIIEIVDALDR_Spectra <- read.csv("I:/EAT_BackUP/Trujillo_34362342_44636349_SN_Data/P2/Complete/Projects/Proteomics/Zymomona/ISPH_3a_FAIMS_infusioon_240K_Is_04iso_infusionSpectra.csv", 
                                           header = TRUE, sep = ",", stringsAsFactors = FALSE)

x <- FAIMS_IspH_AIIEIVDALDR_Spectra

ggplot(x, aes(x=mass.to.charge, y=Intensity)) +
  #geom_line() + 
  #ylim(0,100000)+
  geom_linerange(aes(x=mass.to.charge, ymax=Intensity, ymin=0.75),
                 position = position_jitter(height = 0L, seed = 1L))+
  #scale_y_continuous(breaks = seq(0,60000,10000), limit = c(0,60000), labels = function(x)  format(x, scientific = TRUE)) + #max(x[,3])) +
  #scale_color_manual(values=c('#BB67BC', '#8D7ADD', '#68D6B9', '#54ADD3')) +
  theme_classic()+
  labs(title = "1 min tryptic digest infusion with FAIMS-PRM" ,subtitle = "Strain over-expressinig IspH", y = "Intensity", x ="Time (min)")

##### SPIKED Infusion Spectra extraction of AIEVDQALDR #####

IspH_AIIEIVDALDR_light_heavy <- read.csv("F:/Projects/Proteomics/Zymomona/FAIMS/DataAnalysis/HeavyLight_AIEIVDQALDR_IspHOE_spiked.csv", 
                                                 header = TRUE, sep = ",", stringsAsFactors = FALSE)
colnames(IspH_AIIEIVDALDR_light_heavy) <- c("Time","Intensity","Peptide")

x <- IspH_AIIEIVDALDR_light_heavy
x[,3] <- as.factor(x[,3])
ggplot(x, aes(x=Time, y=Intensity, group = Peptide)) +
  geom_line(aes(color=Peptide), size=2) + 
  #scale_y_continuous(breaks = seq(0,60000,10000), limit = c(0,60000), labels = function(x)  format(x, scientific = TRUE)) + #max(x[,3])) +
  xlim(1.0, 2.0)+
  ylim(0,600000)+
  #scale_color_manual(values=c('#BB67BC', '#8D7ADD', '#68D6B9', '#54ADD3')) +
  theme_classic()+
  labs(title = "1 min tryptic digest infusion with FAIMS-PRM" ,subtitle = "Strain over-expressinig IspH", y = "Intensity", x ="Time (min)")


IspG_ETDIGVTGGGQGK_light_heavy <- read.csv("F:/Projects/Proteomics/Zymomona/FAIMS/DataAnalysis/Infusion_withStandards_Experiment/ETDIGVTGGGQGK_heavy_light_inOEIspG.csv", 
                                         header = TRUE, sep = ",", stringsAsFactors = FALSE)
colnames(IspG_ETDIGVTGGGQGK_light_heavy) <- c("Time","Intensity","Peptide")

x <- IspG_ETDIGVTGGGQGK_light_heavy
x[,3] <- as.factor(x[,3])
ggplot(x, aes(x=Time, y=Intensity, group = Peptide)) +
  geom_line(aes(color=Peptide), size=2) + 
  #scale_y_continuous(breaks = seq(0,60000,10000), limit = c(0,60000), labels = function(x)  format(x, scientific = TRUE)) + #max(x[,3])) +
  xlim(1.0, 2.0)+
  ylim(0,600000)+
  #scale_color_manual(values=c('#BB67BC', '#8D7ADD', '#68D6B9', '#54ADD3')) +
  theme_classic()+
  labs(title = "1 min tryptic digest infusion with FAIMS-PRM" ,subtitle = "Strain over-expressinig IspH", y = "Intensity", x ="Time (min)")



WT_AIIEIVDALDR_light_heavy <- read.csv("F:/Projects/Proteomics/Zymomona/FAIMS/DataAnalysis/HeavyLight_AIEIVDQALDR_WT_spiked.csv", 
                                         header = TRUE, sep = ",", stringsAsFactors = FALSE)
colnames(WT_AIIEIVDALDR_light_heavy) <- c("Time","Intensity","Peptide")

x <- WT_AIIEIVDALDR_light_heavy
x[,3] <- as.factor(x[,3])
ggplot(x, aes(x=Time, y=Intensity, group = Peptide)) +
  geom_line(aes(color=Peptide), size=2) + 
  #scale_y_continuous(breaks = seq(0,60000,10000), limit = c(0,60000), labels = function(x)  format(x, scientific = TRUE)) + #max(x[,3])) +
  xlim(1.0, 2.0)+
  ylim(0,250000)+
  #scale_color_manual(values=c('#BB67BC', '#8D7ADD', '#68D6B9', '#54ADD3')) +
  theme_classic()+
  labs(title = "1 min tryptic digest infusion with FAIMS-PRM" ,subtitle = "Wild Type Strain", y = "Intensity", x ="Time (min)")

##### LC and Spectra extraction of AIEVDQALDR #####

#### Orbi
FAIMS_IspH_AIIEIVDALDR_LCextraction <- read.csv("I:/EAT_BackUP/Trujillo_34362342_44636349_SN_Data/P2/Complete/Projects/Proteomics/Zymomona/20191106_EAT_FAIMS_IspH_1c_extract_AIEIVDQALDR.csv", 
                         header = TRUE, sep = ",", stringsAsFactors = FALSE)
colnames(FAIMS_IspH_AIIEIVDALDR_LCextraction) <- c("Time", "Relative.Abundance")


IspH_MS1_AIIEIVDALDR_Orbi_LCextraction <- read.csv("F:/Projects/Proteomics/Zymomona/FAIMS/DataAnalysis/20200309_ISPH_OrbiMS1_FullChromatogram_45min.csv", 
                                                header = TRUE, sep = ",", stringsAsFactors = FALSE)
colnames(IspH_MS1_AIIEIVDALDR_Orbi_LCextraction) <- c("Time", "Intensity")


IspH_AIIEIVDALDR_Orbi_LCextraction <- read.csv("F:/Projects/Proteomics/Zymomona/FAIMS/DataAnalysis/20200309_ISPH_OrbiMS1_OrbiMS2_FullChromatogram_45min.csv", 
                                               header = TRUE, sep = ",", stringsAsFactors = FALSE)
colnames(IspH_AIIEIVDALDR_Orbi_LCextraction) <- c("Time", "Relative.Abundance")


ggplot(IspH_MS1_AIIEIVDALDR_Orbi_LCextraction, aes(x=Time, y= Intensity)) +
  geom_line()+
  xlim(25.76,26.04)+
  theme_classic()+
  labs(title = "45 min LC - High Resolution MS/MS" ,subtitle = "Retention Window 25.75 min - 26.05 min", y = "Intensity", x ="Retention Time (min)")


#### ION TRAP 
IspH_AIIEIVDALDR_LCextraction <- read.csv("I:/EAT_BackUP/Trujillo_34362342_44636349_SN_Data/P2/Complete/Projects/Proteomics/Zymomona/EAT_Zymo_Isph_3a_20191004_LCextract_AIEIVDQALDR.csv", 
                                                header = TRUE, sep = ",", stringsAsFactors = FALSE)
colnames(IspH_AIIEIVDALDR_LCextraction) <- c("Time", "Relative.Abundance")

ggplot(IspH_AIIEIVDALDR_LCextraction, aes(x=Time, y= Relative.Abundance)) +
  geom_line()+
  xlim(23.23,23.32)+
  theme_classic()+
  labs(title = "45 min LC-MS/MS" ,subtitle = "Retention Window 23.23 min - 23.32 min", y = "Relative Abundance", x ="Retention Time (min)")



FAIMS_IspH_AIIEIVDALDR_SpectraExtraction <- read.csv("I:/EAT_BackUP/Trujillo_34362342_44636349_SN_Data/P2/Complete/Projects/Proteomics/Zymomona/20191106_EAT_FAIMS_IspH_1c_extractSpectra_AIEIVDQALDR.csv", 
                                          header = TRUE, sep = ",", stringsAsFactors = FALSE)
colnames(FAIMS_IspH_AIIEIVDALDR_SpectraExtraction) <- c("mass.to.charge.ratio", "Relative.Abundance")


ggplot(FAIMS_IspH_AIIEIVDALDR_SpectraExtraction, aes(x=mass.to.charge.ratio, y= Relative.Abundance)) +
  #geom_line()+
  geom_linerange(aes(x=mass.to.charge.ratio, ymax=Relative.Abundance, ymin=0.75),
                 position = position_jitter(height = 0L, seed = 1L))+
  theme_classic()+
  labs(title = "MS/MS Spectra" ,subtitle = "XIC of 621.8381 m/z", y = "Relative Abundance", x ="Time (min)")


IspH_AIEIVDALDR_SpectraExtraction <- read.csv("I:/EAT_BackUP/Trujillo_34362342_44636349_SN_Data/P2/Complete/Projects/Proteomics/Zymomona/EAT_Zymo_Isph_3a_20191004_extractSpectra_AIEIVDQALDR.csv", 
                                                     header = TRUE, sep = ",", stringsAsFactors = FALSE)
colnames(IspH_AIEIVDALDR_SpectraExtraction) <- c("mass.to.charge.ratio", "Relative.Abundance")


ggplot(IspH_AIEIVDALDR_SpectraExtraction, aes(x=mass.to.charge.ratio, y= Relative.Abundance)) +
  geom_line()+
  theme_classic()+
  labs(title = "MS/MS Spectra" ,subtitle = "Ion Trap MS/MS scan for precursor 621.8381 m/z", y = "Relative Abundance", x ="Time (min)")

##### OE strain with spiked standard

#Infusion 
ETDIGVTGGGQGK_heavy_light_inOEIspG_infusion <- read.csv("F:/Projects/Proteomics/Zymomona/FAIMS/DataAnalysis/Infusion_withStandards_Experiment/ETDIGVTGGGQGK_heavy_light_inOEIspG.csv", 
                                      header = TRUE, sep = ",", stringsAsFactors = FALSE)

#Spectra _ISPG
ETDIGVTGGGQGK_heavy_light_inOEIspG_spectra <- read.csv("F:/Projects/Proteomics/Zymomona/FAIMS/DataAnalysis/Infusion_withStandards_Experiment/ISPG500ng_250fm_per_uL_heavies_1e6agc_fixmass_LIGHT_HEAVY.csv", 
                                                       header = TRUE, sep = ",", stringsAsFactors = FALSE)
colnames(ETDIGVTGGGQGK_heavy_light_inOEIspG_spectra) <- c("mass.to.charge.ratio", "Intensity", "Ion.Type")
ETDIGVTGGGQGK_heavy_light_inOEIspG_spectra <- ETDIGVTGGGQGK_heavy_light_inOEIspG_spectra[order(ETDIGVTGGGQGK_heavy_light_inOEIspG_spectra$mass.to.charge.ratio),]
ETDIGVTGGGQGK_heavy_light_inOEIspG_spectra$Relative.Abundance <- ETDIGVTGGGQGK_heavy_light_inOEIspG_spectra$Intensity/max(ETDIGVTGGGQGK_heavy_light_inOEIspG_spectra$Intensity)


df <- ETDIGVTGGGQGK_heavy_light_inOEIspG_spectra
background <- df[grep("background",df$Ion.Type),]
transition_light <- df[grep("transition_light",df$Ion.Type),]
transition_heavy <- df[grep("transition_heavy",df$Ion.Type),]

ggplot(background, aes(x=mass.to.charge.ratio, y= Relative.Abundance)) +
  #geom_line()+
  geom_linerange(aes(x=mass.to.charge.ratio, ymax=Relative.Abundance, ymin=0),
                 position = position_jitter(height = 0L, seed = 1L))+
  geom_linerange(data = transition_light, 
                 aes(x=mass.to.charge.ratio, ymax = Relative.Abundance, ymin =0), 
                 position = position_jitter(height = 0L, seed = 1L), color = "red")+
  geom_linerange(data = transition_heavy, 
                 aes(x=mass.to.charge.ratio, ymax = Relative.Abundance, ymin =0), 
                 position = position_jitter(height = 0L, seed = 1L), color = "blue")+
  xlim(150,1250)+
  theme_classic()

ggplot(background, aes(x=mass.to.charge.ratio, y= Relative.Abundance)) +
  #geom_line()+
  geom_linerange(aes(x=mass.to.charge.ratio, ymax=Relative.Abundance, ymin=0),
                 color = "#808080",
                 position = position_jitter(height = 0L, seed = 1L))+
  geom_linerange(data = transition_light, 
                 aes(x=mass.to.charge.ratio, ymax = Relative.Abundance, ymin =0), 
                 position = position_jitter(height = 0L, seed = 1L), color = "red")+
  geom_linerange(data = transition_heavy, 
                 aes(x=mass.to.charge.ratio, ymax = Relative.Abundance, ymin =0), 
                 position = position_jitter(height = 0L, seed = 1L), color = "blue")+
  xlim(550,800)+
  theme_classic()

#Spectra _ISPH
AIEIVDQALDR_heavy_light_inOEIspH_spectra <- read.csv("F:/Projects/Proteomics/Zymomona/FAIMS/DataAnalysis/Infusion_withStandards_Experiment/ISPH500ng_250fm_per_uL_heavies_1e6agc_fixmass_LIGHT_HEAVY.csv", 
                                                       header = TRUE, sep = ",", stringsAsFactors = FALSE)
colnames(AIEIVDQALDR_heavy_light_inOEIspH_spectra) <- c("mass.to.charge.ratio", "Intensity", "Ion.Type")
AIEIVDQALDR_heavy_light_inOEIspH_spectra <- AIEIVDQALDR_heavy_light_inOEIspH_spectra[order(AIEIVDQALDR_heavy_light_inOEIspH_spectra$mass.to.charge.ratio),]
AIEIVDQALDR_heavy_light_inOEIspH_spectra$Relative.Abundance <- AIEIVDQALDR_heavy_light_inOEIspH_spectra$Intensity/max(AIEIVDQALDR_heavy_light_inOEIspH_spectra$Intensity)


df <- AIEIVDQALDR_heavy_light_inOEIspH_spectra
background <- df[grep("background",df$Ion.Type),]
transition_light <- df[grep("transition_light",df$Ion.Type),]
transition_heavy <- df[grep("transition_heavy",df$Ion.Type),]

ggplot(background, aes(x=mass.to.charge.ratio, y= Relative.Abundance)) +
  #geom_line()+
  geom_linerange(aes(x=mass.to.charge.ratio, ymax=Relative.Abundance, ymin=0),
                 position = position_jitter(height = 0L, seed = 1L))+
  geom_linerange(data = transition_light, 
                 aes(x=mass.to.charge.ratio, ymax = Relative.Abundance, ymin =0), 
                 position = position_jitter(height = 0L, seed = 1L), color = "red")+
  geom_linerange(data = transition_heavy, 
                 aes(x=mass.to.charge.ratio, ymax = Relative.Abundance, ymin =0), 
                 position = position_jitter(height = 0L, seed = 1L), color = "blue")+
  xlim(150,1250)+
  theme_classic()

ggplot(background, aes(x=mass.to.charge.ratio, y= Relative.Abundance)) +
  #geom_line()+
  geom_linerange(aes(x=mass.to.charge.ratio, ymax=Relative.Abundance, ymin=0),
                 color = "#808080",
                 position = position_jitter(height = 0L, seed = 1L))+
  geom_linerange(data = transition_light, 
                 aes(x=mass.to.charge.ratio, ymax = Relative.Abundance, ymin =0), 
                 position = position_jitter(height = 0L, seed = 1L), color = "red")+
  geom_linerange(data = transition_heavy, 
                 aes(x=mass.to.charge.ratio, ymax = Relative.Abundance, ymin =0), 
                 position = position_jitter(height = 0L, seed = 1L), color = "blue")+
  xlim(750,1100)+
  theme_classic()




##### Percent TIC explained #####
ProteinProspector_AIEIVQALDR <- read.csv("F:/Projects/Proteomics/Zymomona/FAIMS/DataAnalysis/ProteinProspector_AIEIVDQALDR.csv", 
                                      header = TRUE, sep = ",", stringsAsFactors = FALSE)
ProteinProspector_ETDIGVTGGGQGK <- read.csv("F:/Projects/Proteomics/Zymomona/FAIMS/DataAnalysis/ProteinProspector_ETDIGVTGGGQGK.csv", 
                                         header = TRUE, sep = ",", stringsAsFactors = FALSE)

tic_explained_function <- function(df_proteinProspect, df){
  
  if(length(grep("Low",deparse(substitute(df)))) > 0){
    print(deparse(substitute(df)))
    transitions <- numeric()
    for(i in 1:length(df_proteinProspect$mass.to.charge)){
      transitions <- append(transitions,grep(substr(df_proteinProspect$mass.to.charge[i],1,4), 
                                             substr(df$mass.to.charge.ratio,1,4)))
    }
  }else if(length(grep("High",deparse(substitute(df)))) > 0){
    print(deparse(substitute(df)))
    transitions <- numeric()
    
    for(i in 1:length(df_proteinProspect$mass.to.charge)){
      transitions <- append(transitions,grep(substr(df_proteinProspect$mass.to.charge[i],1,5), 
                                             substr(df$mass.to.charge.ratio,1,5)))
    }
    }else if(length(grep("FAIMS",deparse(substitute(df)))) > 0){
    print(deparse(substitute(df)))
    
    df <- df[-which(df$Relative.Abundance == 0),]
    
    
    transitions <- numeric()
    for(i in 1:length(fragment)){
      transitions <- append(transitions,grep(substr(fragment[i],1,6), 
                                             substr(df$mass.to.charge.ratio,1,6)))
    }
  }
  
  
  df_ions <- df[transitions,]
  
  transitions_reduced <- numeric()
  intensity_toSum <- numeric()
  k=0
  j=0
  for (i in 1:(nrow(df_ions)-1)) {
  #for (i in 1:12) {
    
    if(k>0 & k<=nrow(df_ions)){
      
      
      if(j>=nrow(df_ions)){
        k=i
      }else if(k < nrow(df_ions)){
        k = k+1
        j = k+1
        print(paste0("move to next ion k= ",k, " j= ",j))
      }else if(k == nrow(df_ions)){
        j=k
      }
      
      print(paste0("k=",k, " i=",i))
      print(sub("\\..*", "", df_ions$mass.to.charge.ratio[k]))
      print(sub("\\..*", "", df_ions$mass.to.charge.ratio[j]))
      
      if(sub("\\..*", "", df_ions$mass.to.charge.ratio[k]) == sub("\\..*", "", df_ions$mass.to.charge.ratio[j])){
        abundance_1 <- df_ions$Relative.Abundance[k]
        j = k+1
        abundance_2 <- df_ions$Relative.Abundance[j]
        
        if(abundance_1 < abundance_2){
          print("second ion more intense")
          keep <- df_ions$mass.to.charge.ratio[j]
          keep_intensity <- df_ions$Intensity[j]
          
          
          transitions_reduced <- append(transitions_reduced, keep)
          intensity_toSum <- append(intensity_toSum, keep_intensity)
          
          
          print(transitions_reduced)
          
        }else if(abundance_1 > abundance_2){
          print("first ion more intense")
          
          keep <- df_ions$mass.to.charge.ratio[k]
          keep_intensity <- df_ions$Intensity[k]
          
          if(keep == transitions_reduced[length(transitions_reduced)]){

          }else{
            transitions_reduced <- append(transitions_reduced, keep)
            intensity_toSum <- append(intensity_toSum, keep_intensity)
            }
          
          print(transitions_reduced)
          #k=j+1
          
          if(k == nrow(df_ions)){
            break}else{
              k = k+1
            }
          }
        }else{
          if(j == nrow(df_ions)){
            
            transitions_reduced <- append(transitions_reduced, df_ions$mass.to.charge.ratio[k])
            intensity_toSum <- append(intensity_toSum, df_ions$Intensity[k])
            
            transitions_reduced <- append(transitions_reduced,df_ions$mass.to.charge.ratio[j])
            intensity_toSum <- append(intensity_toSum, df_ions$Intensity[j])
            
            print(transitions_reduced)
            
            if(j >= 27){
              break
            }
          }else{
                        if(transitions_reduced[length(transitions_reduced)] == df_ions$mass.to.charge.ratio[k]){
              j = k+1
              transitions_reduced <- append(transitions_reduced,df_ions$mass.to.charge.ratio[j])
              intensity_toSum <- append(intensity_toSum, df_ions$Intensity[j])
              
            }else{
            transitions_reduced <- append(transitions_reduced,df_ions$mass.to.charge.ratio[k])
            intensity_toSum <- append(intensity_toSum, df_ions$Intensity[k])
            }
          }
            print(transitions_reduced)
        }
      

      }else if(k == 0){
      print(paste0("i=",i))
      
      #print(sub("\\..*", "", df_ions$mass.to.charge.ratio[i]))
      #print(sub("\\..*", "", df_ions$mass.to.charge.ratio[i+1]))
      
      if(sub("\\..*", "", df_ions$mass.to.charge.ratio[i]) == sub("\\..*", "", df_ions$mass.to.charge.ratio[i+1])){
        abundance_1 <- df_ions$Relative.Abundance[i]
        j = i+1
        abundance_2 <- df_ions$Relative.Abundance[j]
        
        if(abundance_1 < abundance_2){
          print("second ion more intense")
          keep <- df_ions$mass.to.charge.ratio[j]
          keep_intensity <- df_ions$Intensity[j]
          transitions_reduced <- append(transitions_reduced, keep)
          intensity_toSum <- append(intensity_toSum, keep_intensity)
          print(transitions_reduced)
          
          if(sub("\\..*", "", df_ions$mass.to.charge.ratio[j]) == sub("\\..*", "", df_ions$mass.to.charge.ratio[j+1])){
            k=0
          }else{
            k=j
          }
          
        }else if(abundance_1 > abundance_2){
          print("first ion more intense when K is zero")
          keep <- df_ions$mass.to.charge.ratio[i]
          keep_intensity <- df_ions$Intensity[i]
          
          if(i > 1){
            if(keep == transitions_reduced[i-1]){
              
            }else{
              transitions_reduced <- append(transitions_reduced, keep)
              intensity_toSum <- append(intensity_toSum, keep_intensity)
              }
            
          }else{
            transitions_reduced <- append(transitions_reduced, keep)
            intensity_toSum <- append(intensity_toSum, keep_intensity)
          }
          print(transitions_reduced)
        
          k=j
          
        }else if (abundance_1 == abundance_2){
          
          k=j
        }
      }else{ 
        if(k==0){
          transitions_reduced <- append(transitions_reduced,df_ions$mass.to.charge.ratio[i])
          intensity_toSum <- append(intensity_toSum, df_ions$Intensity[i])
          print(transitions_reduced)
        }else{
          transitions_reduced <- append(transitions_reduced,df_ions$mass.to.charge.ratio[k])
          intensity_toSum <- append(intensity_toSum, df_ions$Intensity[k])
          k=0
        }
        
      }
      }

    }
  
    
  Tic_explained <- (sum(intensity_toSum[-c(16,17,18)]))/(sum(df$Intensity))*100
  
  return(Tic_explained)
  
  }



# df_proteinProspect <- ProteinProspector_AIEIVQALDR
# df <- noFAIMS_AIIEIVDALDR_Spectra_RelativeAbundance
Tic_explained_High_Res_AIEIDVQALDR <- tic_explained_function(ProteinProspector_AIEIVQALDR, High_Res_IspH_AIIEIVDALDR_Spectra.RelativeAbundance)

Tic_explained_Low_Res_AIEIDVQALDR <- tic_explained_function(ProteinProspector_AIEIVQALDR,Low_Res_IspH_AIIEIVDALDR_Spectra.RelaiveAbundance)

Tic_explained_High_Res_ETDIGVTGGGQGK <- tic_explained_function(ProteinProspector_ETDIGVTGGGQGK, High_Res_IspG_ETDIGVTGGGGQGK_Spectra.RelativeAbundance)

Tic_explained_Low_Res_ETDIGVTGGGQGK <- tic_explained_function(ProteinProspector_ETDIGVTGGGQGK, Low_Res_IspG_ETDIGVTGGGGQGK_Spectra.RelativeAbundance)

Tic_explained_FAIMS_AIEIDVQALDR <- tic_explained_function(ProteinProspector_AIEIVQALDR,FAIMS_AIIEIVDALDR_Spectra_RelativeAbundance)

Tic_explained_noFAIMS_AIEIDVQALDR <- tic_explained_function(ProteinProspector_AIEIVQALDR, noFAIMS_AIIEIVDALDR_Spectra_RelativeAbundance)

# Peptide abundances
ISPH_Jan25_1 <- read.csv("I:/EAT_BackUP/Trujillo_34362342_44636349_SN_Data/P2/Complete/Projects/Proteomics/Zymomona/DataAnalysis/SkyLine_LibraryComparison/ISPH/IW1_240K_MI1_AGC1e06_noFAIMS_1_20200125.csv", 
                            header = TRUE, sep = ",", stringsAsFactors = FALSE)
ISPH_Jan25_2 <- read.csv("I:/EAT_BackUP/Trujillo_34362342_44636349_SN_Data/P2/Complete/Projects/Proteomics/Zymomona/DataAnalysis/SkyLine_LibraryComparison/ISPH/IW1_240K_MI1_AGC1e06_noFAIMS_2_20200125.csv", 
                         header = TRUE, sep = ",", stringsAsFactors = FALSE)
ISPH_Jan25_3 <- read.csv("I:/EAT_BackUP/Trujillo_34362342_44636349_SN_Data/P2/Complete/Projects/Proteomics/Zymomona/DataAnalysis/SkyLine_LibraryComparison/ISPH/IW1_240K_MI1_AGC1e06_noFAIMS_3_20200125.csv", 
                         header = TRUE, sep = ",", stringsAsFactors = FALSE)

ISPH_Jan25_1to3 <- read.csv("I:/EAT_BackUP/Trujillo_34362342_44636349_SN_Data/P2/Complete/Projects/Proteomics/Zymomona/DataAnalysis/SkyLine_LibraryComparison/ISPH/IW1_240K_MI1_AGC1e06_noFAIMS_1-3_20200125.csv", 
                         header = TRUE, sep = ",", stringsAsFactors = FALSE)

ISPG_Jan25_1to3 <- read.csv("I:/EAT_BackUP/Trujillo_34362342_44636349_SN_Data/P2/Complete/Projects/Proteomics/Zymomona/DataAnalysis/SkyLine_LibraryComparison/ISPG/IW1_240K_MI1_AGC1e06_noFAIMS_1-3_ISPGPeptides_20200125.csv", 
                            header = TRUE, sep = ",", stringsAsFactors = FALSE)

ISPH_plot(ISPH_Jan25_1)
ISPH_plot(ISPH_Jan25_2)
ISPH_plot(ISPH_Jan25_3)
ISPH_plot(ISPH_Jan25_1to3)
ISPG_plot(ISPG_Jan25_1to3)

# NO-FAIMS Fragmentation coverage
AIEIVDQALDR_ISPH <- read.csv("I:/EAT_BackUP/Trujillo_34362342_44636349_SN_Data/P2/Complete/Projects/Proteomics/Zymomona/DataAnalysis/SkyLine_LibraryComparison/ISPH/AIEIVDQALDR_libraryComparison.csv", 
                         header = TRUE, sep = ",", stringsAsFactors = FALSE)
GLPVVDATCPLVNK_ISPH <- read.csv("I:/EAT_BackUP/Trujillo_34362342_44636349_SN_Data/P2/Complete/Projects/Proteomics/Zymomona/DataAnalysis/SkyLine_LibraryComparison/ISPH/GLPVVDATCPLVNK_libraryComparison.csv", 
                                header = TRUE, sep = ",", stringsAsFactors = FALSE)
FTDVIGPDTSDICYATQNR_ISPH <- read.csv("I:/EAT_BackUP/Trujillo_34362342_44636349_SN_Data/P2/Complete/Projects/Proteomics/Zymomona/DataAnalysis/SkyLine_LibraryComparison/ISPH/FTDVIGPDTSDICYATQNR_libraryComparison.csv", 
                                     header = TRUE, sep = ",", stringsAsFactors = FALSE)

VSLSADPEQEVR_ISPG <- read.csv("I:/EAT_BackUP/Trujillo_34362342_44636349_SN_Data/P2/Complete/Projects/Proteomics/Zymomona/DataAnalysis/SkyLine_LibraryComparison/ISPG/VSLSADPEQEVR_libraryComparison.csv", 
                              header = TRUE, sep = ",", stringsAsFactors = FALSE)

ETDIGVTGGGQGK_ISPG <- read.csv("I:/EAT_BackUP/Trujillo_34362342_44636349_SN_Data/P2/Complete/Projects/Proteomics/Zymomona/DataAnalysis/SkyLine_LibraryComparison/ISPG/ETDIGVTGGGQGK_libraryComparison.csv", 
                               header = TRUE, sep = ",", stringsAsFactors = FALSE)

AVDCPLHLGITEAGGLIGGTVK_ISPG <- read.csv("I:/EAT_BackUP/Trujillo_34362342_44636349_SN_Data/P2/Complete/Projects/Proteomics/Zymomona/DataAnalysis/SkyLine_LibraryComparison/ISPG/AVDCPLHLGITEAGGLIGGTVK_libraryComparison.csv", 
                                        header = TRUE, sep = ",", stringsAsFactors = FALSE)

# FAIMS Fragmentation coverage
FAIMS_AIEIVDQALDR_ISPH <- read.csv("I:/EAT_BackUP/Trujillo_34362342_44636349_SN_Data/P2/Complete/Projects/Proteomics/Zymomona/DataAnalysis/SkyLine_LibraryComparison/ISPH/FAIMS_AIEIVDQALDR_libraryComparison.csv", 
                             header = TRUE, sep = ",", stringsAsFactors = FALSE)
FAIMS_GLPVVDATCPLVNK_ISPH <- read.csv("I:/EAT_BackUP/Trujillo_34362342_44636349_SN_Data/P2/Complete/Projects/Proteomics/Zymomona/DataAnalysis/SkyLine_LibraryComparison/ISPH/FAIMS_GLPVVDATCPLVNK_libraryComparison.csv", 
                                header = TRUE, sep = ",", stringsAsFactors = FALSE)
FAIMS_FTDVIGPDTSDICYATQNR_ISPH <- read.csv("I:/EAT_BackUP/Trujillo_34362342_44636349_SN_Data/P2/Complete/Projects/Proteomics/Zymomona/DataAnalysis/SkyLine_LibraryComparison/ISPH/FAIMS_FTDVIGPDTSDICYATQNR_libraryComparison.csv", 
                                     header = TRUE, sep = ",", stringsAsFactors = FALSE)

FAIMS_VSLSADPEQEVR_ISPG <- read.csv("I:/EAT_BackUP/Trujillo_34362342_44636349_SN_Data/P2/Complete/Projects/Proteomics/Zymomona/DataAnalysis/SkyLine_LibraryComparison/ISPG/FAIMS_VSLSADPEQEVR_libraryComparison.csv", 
                              header = TRUE, sep = ",", stringsAsFactors = FALSE)

FAIMS_ETDIGVTGGGQGK_ISPG <- read.csv("I:/EAT_BackUP/Trujillo_34362342_44636349_SN_Data/P2/Complete/Projects/Proteomics/Zymomona/DataAnalysis/SkyLine_LibraryComparison/ISPG/FAIMS_ETDIGVTGGGQGK_libraryComparison.csv", 
                               header = TRUE, sep = ",", stringsAsFactors = FALSE)



## For pheatmap_1.0.8 and later:rotate column names by 0 degrees
draw_colnames_90 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 0, gp = gpar(...))
  return(res)}

peptide_coverage <- function(x){
  m <- as.matrix(x)
  
  dimnames(m) <- list(c("PROSIT","HighRes", "LowRes"), names(x)) 
  show(m)
  
  
  ## 'Overwrite' default draw_colnames with your own version 
  assignInNamespace(x="draw_colnames", value="draw_colnames_90",
                     ns=asNamespace("pheatmap"))
  
  h <- pheatmap(m,
                breaks = c(0,0.9,2),
                color = c( "white", "#A4B8C4"),
                border_color = "gray",
                cluster_cols = FALSE,
                cluster_rows = FALSE)
  

  return(h)
  
}

# Size of image is 536 X 100

#AIEIVDQALDR_Prosit_HighRes_LowRes.svg
AIEIVDQALDR <- peptide_coverage(AIEIVDQALDR_ISPH)

#GLPVVDATCPLVNK_Prosit_HighRes_LowRes.svg
GLPVVDATCPLVNK <- peptide_coverage(GLPVVDATCPLVNK_ISPH)

#FTDVIGPDTSDICYATQNR_Prosit_HighRes_LowRes.svg
FTDVIGPDTSDICYATQNR <- peptide_coverage(FTDVIGPDTSDICYATQNR_ISPH)

#VSLSADPEQEVR_Prosit_HighRes_LowRes.svg
VSLSADPEQEVR <- peptide_coverage(VSLSADPEQEVR_ISPG)

#ETDIGVTGGGQGK_Prosit_HighRes_LowRes.svg
ETDIGVTGGGQGK <- peptide_coverage(ETDIGVTGGGQGK_ISPG)

#AVDCPLHLGITEAGGLIGGTVK_Prosit_HighRes_LowRes.svg
AVDCPLHLGITEAGGLIGGTVK <- peptide_coverage(AVDCPLHLGITEAGGLIGGTVK_ISPG)

### FAIMS

#FAIMS_AIEIVDQALDR_Prosit_HighRes_LowRes.svg
FAIMS_AIEIVDQALDR <- peptide_coverage(FAIMS_AIEIVDQALDR_ISPH)

#FAIMS_ETDIGVTGGGQGK_Prosit_HighRes_LowRes.svg
FAIMS_ETDIGVTGGGQGK <- peptide_coverage(FAIMS_ETDIGVTGGGQGK_ISPG)


##### Standards ####

Standards_IspG <- read.csv("F:/Projects/Proteomics/Zymomona/FAIMS/DataAnalysis/PeakArea_Standard_IspG.csv", 
                             header = TRUE, sep = ",", stringsAsFactors = FALSE)

Standards_IspH <- read.csv("F:/Projects/Proteomics/Zymomona/FAIMS/DataAnalysis/PeakArea_Standard_IspH_FAIMS.csv", 
                           header = TRUE, sep = ",", stringsAsFactors = FALSE)

ggplot(Standards_IspG, aes(x=log2(Concentration), y= y9._768.4090)) +
  geom_point(size = 3) +
  geom_line(linetype = 2)+
  geom_point(data = Standards_IspG, aes(x=log2(Concentration), y=y7_612.3191), size = 3, color = "red")+
  geom_line(data = Standards_IspG, aes(x=log2(Concentration), y=y7_612.3191), linetype = 2, color = "red")+
  geom_point(data = Standards_IspG, aes(x=log2(Concentration), y=y6_511.2714), size = 3, color = "green")+
  geom_line(data = Standards_IspG, aes(x=log2(Concentration), y=y6_511.2714), linetype = 2, color = "green")+
  theme_light()+
  labs(title = "Standard IspG" ,subtitle = "Transitions y6, y7, and y9", y = "Log2(Intensity)", x ="Log2(Concentration - attomole)")
  

ggplot(Standards_IspH, aes(x=log2(Concentration), y=y6_727.3609)) +
  geom_point(size = 3) +
  geom_line(linetype = 2)+
  geom_point(data = Standards_IspH, aes(x=log2(Concentration), y=y7_826.4293), size = 3, color = "red")+
  geom_line(data = Standards_IspH, aes(x=log2(Concentration), y=y7_826.4293), linetype = 2, color = "red")+
  geom_point(data = Standards_IspH, aes(x=log2(Concentration), y=y9_1068.5559), size = 3, color = "green")+
  geom_line(data = Standards_IspH, aes(x=log2(Concentration), y=y9_1068.5559), linetype = 2, color = "green")+
  theme_light()+
    labs(title = "Standard IspH - FAIMS" ,subtitle = "Transitions y6, y7,and y9", y = "Log2(Intensity)", x ="Log2(Concentration) (femtomole)") 
  

#### Infusion Experiment ####
#FAIMS
FAIMS_Infusion_Experiment <- read.csv("F:/Projects/Proteomics/Zymomona/FAIMS/DataAnalysis/Ratio_LightToHeavy_FAIMS.csv", 
                                        header = TRUE, sep = ",", stringsAsFactors = FALSE)

meta <- FAIMS_Infusion_Experiment
meta$Sample <- gsub("-", "_", FAIMS_Infusion_Experiment$Sample)
meta$Sample <- sub(" Ratio.*", "", meta$Sample)

meta$Sample <- ifelse(substr(meta$Sample, 1, 1) == "_", sub("^.", "", meta$Sample), "")
meta_sample_exanded <- read.table(text = meta$Sample, sep = "_", colClasses = "character")
colnames(meta_sample_exanded) <- c("Strain", "Replicate", "Time", "Index")
meta_sample_exanded$Time <- as.numeric(sub("7point5", "7.5", meta_sample_exanded$Time))

meta_sample_exanded$TimeReplicates <- paste0(meta_sample_exanded$Strain, "_", meta_sample_exanded$Time)


FAIMS_Infusion_Experiment <- cbind(meta_sample_exanded$TimeReplicates, FAIMS_Infusion_Experiment)

average_time.point <- aggregate(FAIMS_Infusion_Experiment[,3:ncol(FAIMS_Infusion_Experiment)],
                                list(FAIMS_Infusion_Experiment$`meta_sample_exanded$TimeReplicates`), mean)

std_time.point_extended <- aggregate(FAIMS_Infusion_Experiment[,3:ncol(FAIMS_Infusion_Experiment)],
                                list(FAIMS_Infusion_Experiment$`meta_sample_exanded$TimeReplicates`), sd)

# order
average_time.point <- unique(average_time.point[match(meta_sample_exanded$TimeReplicates, average_time.point$Group.1),])
colnames(average_time.point) <- c("TimeReplicates", "AIEIVDQALDR", "ETDIGVTGGGQGK")

std_time.point <- unique(std_time.point_extended[match(meta_sample_exanded$TimeReplicates, std_time.point_extended$Group.1),])
colnames(std_time.point) <- c("TimeReplicates", "AIEIVDQALDR", "ETDIGVTGGGQGK")

meta_average_time.point <- as.data.frame(stringr::str_split_fixed(average_time.point$TimeReplicates,"_",2))
average_time.point_extended <- cbind(meta_average_time.point, average_time.point)
average_time.point_extended$V2 <- as.numeric(as.character(average_time.point_extended$V2))
colnames(average_time.point_extended) <- c("Strain", "Time", "TimeReplicates", "AIEIVDQALDR", "ETDIGVTGGGQGK")

std_time.point_extended <- cbind(meta_average_time.point, std_time.point)
std_time.point_extended$V2 <- as.numeric(as.character(std_time.point_extended$V2))
colnames(std_time.point_extended) <- c("Strain", "Time", "TimeReplicates", "AIEIVDQALDR", "ETDIGVTGGGQGK")

avg_std_time.point_extended <- cbind(average_time.point_extended,std_time.point_extended[,4:5])
colnames(avg_std_time.point_extended) <- c("Strain", "Time", "TimeReplicates", "AIEIVDQALDR", "ETDIGVTGGGQGK", "error.AIEIVDQALDR", "error.ETDIGVTGGGQGK")
# 
# average_time.point_extended <- average_time.point_extended %>%
#   mutate(Time = factor(Time, level = seq(0,180,1))) 

#### linear model ####

linearMod <- lm(AIEIVDQALDR ~ Strain + Time + Strain*Time, data = avg_std_time.point_extended)
plot(linearMod)
anova1 <- aov(AIEIVDQALDR ~ Strain + Time + Strain*Time, data = avg_std_time.point_extended)
summary(anova1)
anova1_lm <- aov(linearMod)

#IspH
linearMod_strain_time_ETDIGVTGGGQGK <- lm(ETDIGVTGGGQGK ~ Strain + Time , data = avg_std_time.point_extended)
plot(linearMod_strain_time_ETDIGVTGGGQGK)
anova_strain_time_ETDIGVTGGGQGK <- aov(linearMod_strain_time_ETDIGVTGGGQGK)
summary(anova_strain_time_ETDIGVTGGGQGK)

TukeyHSD(anova_strain_time_ETDIGVTGGGQGK)

#IspG
linearMod_strain_time <- lm(AIEIVDQALDR ~ Strain + Time , data = avg_std_time.point_extended)
plot(linearMod_strain_time)
anova_strain_time <- aov(linearMod_strain_time)
summary(anova_strain_time)

TukeyHSD(anova_strain_time)

average_time.point_extended$Strain <- as.factor(average_time.point_extended$Strain)

ggplot(average_time.point_extended, aes(x= Time, y=AIEIVDQALDR, group=Strain)) +
  geom_point(aes(color = Strain)) +
  geom_line(aes(color = Strain)) +
  xlim(0,180)+
  theme_classic()+
  labs(title = "Time Course IspH Over Expression Monitoring", subtitle = "Light to Heavy Ratio - AIEIVDQALDR",
       x = "Time (min)", y = "Light to Heavy Ratio")


ggplot(avg_std_time.point_extended, aes(x= Time, y=AIEIVDQALDR, group=Strain)) +
  geom_point(aes(color = Strain)) +
  geom_line(aes(color = Strain)) +
  geom_errorbar(aes(ymin=AIEIVDQALDR-error.AIEIVDQALDR, ymax=AIEIVDQALDR+error.AIEIVDQALDR),
                alpha = 0.5,
                width = 0,
                size= 10,
                color= "#6D696F") +
  xlim(0,180)+
  theme_classic()+
  labs(title = "Time Course IspH Over Expression Monitoring", subtitle = "Light to Heavy Ratio - AIEIVDQALDR",
       x = "Time (min)", y = "Light to Heavy Ratio")

ggplot(avg_std_time.point_extended, aes(x= Time, y=AIEIVDQALDR, group=Strain)) +
  geom_point(aes(color = Strain)) +
  geom_line(aes(color = Strain)) +
  geom_ribbon(aes(ymin = AIEIVDQALDR-error.AIEIVDQALDR, ymax = AIEIVDQALDR+error.AIEIVDQALDR), 
              alpha = 0.1,
              linetype =1,
              size =1)+
  xlim(0,180)+
  theme_classic()+
  labs(title = "Time Course IspH Over Expression Monitoring", subtitle = "Light to Heavy Ratio - AIEIVDQALDR",
       x = "Time (min)", y = "Light to Heavy Ratio")

ggplot(average_time.point_extended[1:7,], aes(x= Time, y=AIEIVDQALDR)) +
  geom_line() +
  stat_smooth()+
  xlim(0,180)+
  theme_classic()+
  labs(title = "Time Course IspH Over Expression Monitoring", subtitle = "Light to Heavy Ratio - AIEIVDQALDR",
       x = "Time (min)", y = "Light to Heavy Ratio")

average_time.point_extended$Strain <- as.factor(average_time.point_extended$Strain)
ggplot(average_time.point_extended, aes(x= Time, y=ETDIGVTGGGQGK, group=Strain)) +
  geom_point(aes(color = Strain)) +
  geom_line(aes(color = Strain)) +
  theme_classic()

ggplot(avg_std_time.point_extended, aes(x= Time, y=ETDIGVTGGGQGK, group=Strain)) +
  geom_point(aes(color = Strain)) +
  geom_line(aes(color = Strain)) +
  geom_errorbar(aes(ymin=ETDIGVTGGGQGK-error.ETDIGVTGGGQGK, ymax=ETDIGVTGGGQGK+error.ETDIGVTGGGQGK),
                alpha = 0.1,
                width = 0,
                size= 10,
                color= "#6D696F") +
  xlim(0,180)+
  theme_classic()+
  labs(title = "Time Course IspG Over Expression Monitoring", subtitle = "Light to Heavy Ratio - ETDIGVTGGGQGK",
       x = "Time (min)", y = "Light to Heavy Ratio")
