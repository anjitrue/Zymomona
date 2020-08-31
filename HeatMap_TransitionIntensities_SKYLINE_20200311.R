library(pheatmap)
library(RColorBrewer)

#### Upload Data ####

palette(c("#5ABCB2", "#D6F8D6","#D6F8D6" ,"mediumorchid2","darkgoldenrod1", 
          "hotpink3", "red2", "sienna2","slategray4","mediumturquoise", 
          "deepskyblue", "orangered", "midnightblue"))

# Read in SkyLine output proteingroups, I filtered out the contaminants and the fasta header overspill

# data.frame of Eclipse 45 min runs using HP column. Searched with new aliqouts to confirm the strange intensity trends 
VSLSADPEQEVR_libraryComparison <- read.csv("H:/Projects/Proteomics/Zymomona/FAIMS/Zymo/LC_runs/VSLSADPEQEVR_libraryComparison.csv", 
                                                               header = TRUE, sep = ",", stringsAsFactors = FALSE)
rownames(VSLSADPEQEVR_libraryComparison) <- VSLSADPEQEVR_libraryComparison[,1]
VSLSADPEQEVR_libraryComparison <- VSLSADPEQEVR_libraryComparison[,-1]


cols <- brewer.pal(3, "BuGn")
scaleRYG <- colorRampPalette(cols)(11)

pheatmap(VSLSADPEQEVR_libraryComparison,
         color = scaleRYG,
         cluster_cols = FALSE,
         scale = "column")

                                                            
                                                            