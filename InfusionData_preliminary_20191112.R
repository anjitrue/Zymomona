library(ggplot2)
library(plotly)
library(pheatmap)

#### WT_3c_FAIMS_ISPHpeptides_infusion_20191112 ####

WT_FAIMS_ISPH <- read.csv("G:/Projects/Proteomics/Zymomona/WT_3c_FAIMS_ISPHpeptides_infusion_20191112.csv", 
                                            header = TRUE, sep = ",", stringsAsFactors = FALSE)

WT_FAIMS_ISPH$Peptide <- as.factor(WT_FAIMS_ISPH$Peptide)

p <- ggplot(WT_FAIMS_ISPH, aes(x=Retention.Time, y=Intensity, group = Peptide)) +
  geom_line(aes(color=Peptide), size=2) + 
  scale_color_manual(values=c('#BB67BC', '#8D7ADD', '#68D6B9', '#54ADD3'))  
p




#### ISPH_3a_FAIMS_ISPHpeptides_infusion_20191112 ####

ISPH_FAIMS_ISPH <- read.csv("G:/Projects/Proteomics/Zymomona/ISPH_3a_FAIMS_ISPHpeptides_infusion_20191112.csv", 
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

#### ISPH_3a_FAIMS_ISPHpeptides_infusion_240K ####

ISPH_FAIMS_240K <- read.csv("G:/Projects/Proteomics/Zymomona/ISPH_3a_FAIMS_ISPHpeptides_infusion_240K.csv", 
                            header = TRUE, sep = ",", stringsAsFactors = FALSE)

ISPH_plot(ISPH_FAIMS_240K)

##### Library Comparison ####

# Peptide abundances
ISPH_Jan25_1 <- read.csv("G:/Projects/Proteomics/Zymomona/DataAnalysis/SkyLine_LibraryComparison/ISPH/IW1_240K_MI1_AGC1e06_noFAIMS_1_20200125.csv", 
                            header = TRUE, sep = ",", stringsAsFactors = FALSE)
ISPH_Jan25_2 <- read.csv("G:/Projects/Proteomics/Zymomona/DataAnalysis/SkyLine_LibraryComparison/ISPH/IW1_240K_MI1_AGC1e06_noFAIMS_2_20200125.csv", 
                         header = TRUE, sep = ",", stringsAsFactors = FALSE)
ISPH_Jan25_3 <- read.csv("G:/Projects/Proteomics/Zymomona/DataAnalysis/SkyLine_LibraryComparison/ISPH/IW1_240K_MI1_AGC1e06_noFAIMS_3_20200125.csv", 
                         header = TRUE, sep = ",", stringsAsFactors = FALSE)

ISPH_Jan25_1to3 <- read.csv("G:/Projects/Proteomics/Zymomona/DataAnalysis/SkyLine_LibraryComparison/ISPH/IW1_240K_MI1_AGC1e06_noFAIMS_1-3_20200125.csv", 
                         header = TRUE, sep = ",", stringsAsFactors = FALSE)

ISPG_Jan25_1to3 <- read.csv("G:/Projects/Proteomics/Zymomona/DataAnalysis/SkyLine_LibraryComparison/ISPG/IW1_240K_MI1_AGC1e06_noFAIMS_1-3_ISPGPeptides_20200125.csv", 
                            header = TRUE, sep = ",", stringsAsFactors = FALSE)

ISPH_plot(ISPH_Jan25_1)
ISPH_plot(ISPH_Jan25_2)
ISPH_plot(ISPH_Jan25_3)
ISPH_plot(ISPH_Jan25_1to3)
ISPG_plot(ISPG_Jan25_1to3)

# Fragmentation coverage
AIEIVDQALDR_ISPH <- read.csv("G:/Projects/Proteomics/Zymomona/DataAnalysis/SkyLine_LibraryComparison/ISPH/AIEIVDQALDR_libraryComparison.csv", 
                         header = TRUE, sep = ",", stringsAsFactors = FALSE)
GLPVVDATCPLVNK_ISPH <- read.csv("G:/Projects/Proteomics/Zymomona/DataAnalysis/SkyLine_LibraryComparison/ISPH/GLPVVDATCPLVNK_libraryComparison.csv", 
                                header = TRUE, sep = ",", stringsAsFactors = FALSE)
FTDVIGPDTSDICYATQNR_ISPH <- read.csv("G:/Projects/Proteomics/Zymomona/DataAnalysis/SkyLine_LibraryComparison/ISPH/FTDVIGPDTSDICYATQNR_libraryComparison.csv", 
                                     header = TRUE, sep = ",", stringsAsFactors = FALSE)

VSLSADPEQEVR_ISPG <- read.csv("G:/Projects/Proteomics/Zymomona/DataAnalysis/SkyLine_LibraryComparison/ISPG/VSLSADPEQEVR_libraryComparison.csv", 
                              header = TRUE, sep = ",", stringsAsFactors = FALSE)

ETDIGVTGGGQGK_ISPG <- read.csv("G:/Projects/Proteomics/Zymomona/DataAnalysis/SkyLine_LibraryComparison/ISPG/ETDIGVTGGGQGK_libraryComparison.csv", 
                               header = TRUE, sep = ",", stringsAsFactors = FALSE)

AVDCPLHLGITEAGGLIGGTVK_ISPG <- read.csv("G:/Projects/Proteomics/Zymomona/DataAnalysis/SkyLine_LibraryComparison/ISPG/AVDCPLHLGITEAGGLIGGTVK_libraryComparison.csv", 
                                        header = TRUE, sep = ",", stringsAsFactors = FALSE)

peptide_coverage <- function(x){
  m <- as.matrix(x)
  
  dimnames(m) <- list(c("HighRes", "LowRes"), names(x)) 
  show(m)
  
  h <- pheatmap(m,
                breaks = c(0,0.9,2),
                color = c( "white", "#A4B8C4"),
                border_color = "gray",
                cluster_cols = FALSE,
                cluster_rows = FALSE)

  return(h)
  
}

AIEIVDQALDR <- peptide_coverage(AIEIVDQALDR_ISPH)

GLPVVDATCPLVNK <- peptide_coverage(GLPVVDATCPLVNK_ISPH)

FTDVIGPDTSDICYATQNR <- peptide_coverage(FTDVIGPDTSDICYATQNR_ISPH)

VSLSADPEQEVR <- peptide_coverage(VSLSADPEQEVR_ISPG)

ETDIGVTGGGQGK <- peptide_coverage(ETDIGVTGGGQGK_ISPG)

AVDCPLHLGITEAGGLIGGTVK <- peptide_coverage(AVDCPLHLGITEAGGLIGGTVK_ISPG)

