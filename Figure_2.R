

#'[-----------------------------------------------------------]
#'[-----------------------------------------------------------]
#'[ Mass-spec dataset analysis ]


library("dplyr")
library("DEP")
library("tidyr")
library("purrr")
library("ggplot2")
library("SummarizedExperiment")


getwd()
#'[-----Here we converted input log2 intenisty into base 0 original state-------------------------------------------------------]
data_unique <- read.csv("Nes_vs_zuc_tom20__fixed_input_20210723_with_log2intensity_to_base0.csv", header = TRUE)
data_unique[1:2, 12:24]
#Organism Length       name         ID Intensity.NES_TurboID_1 Intensity.NES_TurboID_2
#1 Drosophila melanogaster (Fruit fly)    276   l(2)37Cc A0A023GQA5               165028020               158184801
#2 Drosophila melanogaster (Fruit fly)    147 A0A023T5E7 A0A023T5E7                      NA                      NA
#Intensity.NES_TurboID_3 Intensity.Tom20_TurboID_1 Intensity.Tom20_TurboID_2 Intensity.Tom20_TurboID_3
#1               159202139                6467619324                6496327549                8035495985
#2                      NA                 461927918                 512452896                 450279745
#Intensity.Zuc_TurboID_1 Intensity.Zuc_TurboID_2 Intensity.Zuc_TurboID_3
#1              1807521498              1255997764              1669675885
#2               220582849               181551355               173827743

dim(data_unique) #1866   21


#'[------------------------------------------------------------]
getwd() 

UbiLength_ExpDesign <- read.csv("summarize_object_deg_mass_spec_zuc_nes_20210719.csv", header = TRUE)
UbiLength_ExpDesign
#            label condition replicate
#1   NES_TurboID_1      Ctrl         1
#2   NES_TurboID_2      Ctrl         2
#3   NES_TurboID_3      Ctrl         3
#4 Tom20_TurboID_1     Tom20         1
#5 Tom20_TurboID_2     Tom20         2
#6 Tom20_TurboID_3     Tom20         3
#7   Zuc_TurboID_1       Zuc         1
#8   Zuc_TurboID_2       Zuc         2
#9   Zuc_TurboID_3       Zuc         3


LFQ_columns <- grep("TurboID", colnames(data_unique)) # get LFQ column numbers
experimental_design <- UbiLength_ExpDesign
data_se <- make_se(data_unique, LFQ_columns, experimental_design)


# Generate a SummarizedExperiment object by parsing condition information from the column names
LFQ_columns <- grep("TurboID", colnames(data_unique)) # get LFQ column numbers
data_se_parsed <- make_se_parse(data_unique, LFQ_columns)

# Let's have a look at the SummarizedExperiment object
data_se

#class: SummarizedExperiment
#dim: 1860 9
#metadata(0):
#  assays(1): ''
#rownames(1860): zuc sgg ... Smn 0.108
#rowData names(6): Protein Proteinnames ... name ID
#colnames(9): Ctrl_1 Ctrl_2 ... Zuc_2 Zuc_3
#colData names(4): label ID condition replicate

#'[Filter on missing values]

# Plot a barplot of the protein identification overlap between samples
plot_frequency(data_se)

#'[This leaves our dataset with missing values, which need to be imputed.]

# Filter for proteins that are identified in all replicates of at least one condition
data_filt <- filter_missval(data_se, thr = 0)

# After filtering, the number of identified proteins per sample can be plotted as well as the overlap in identifications between samples

# Plot a barplot of the number of identified proteins per samples
plot_numbers(data_filt)

# Plot a barplot of the protein identification overlap between samples
plot_coverage(data_filt)

#'[Normalization]
#'[vsn]

# Normalize the data
data_norm <- normalize_vsn(data_filt)

# Visualize normalization by boxplots for all samples before and after normalization
plot_normalization(data_filt, data_norm)

#'[Impute data for missing values]

# Plot a heatmap of proteins with missing values
plot_missval(data_filt)

# Plot intensity distributions and cumulative fraction of proteins with and without missing values
plot_detect(data_filt)

#'[Here we only focused on minDET according to prof request!]
data_impMinDet <- impute(data_norm, fun = "MinDet", q = 0.01)


?impute()


#'[The effect of the imputation on the distributions can be visualized.]

plot_imputation(data_norm, data_impMinDet)

#'[Differential enrichment analysis]

#'[Here we only focused on minDET according to prof request!]
#'[MinDET germ cell]
data_diff_manual <- test_diff(data_impMinDet, type = "manual",
                              test = c("Tom20_vs_control", "Zuc_vs_control", 'Zuc_vs_Tom20'))

#'[Tested contrasts: Tom20_vs_control, Zuc_vs_control, Zuc_vs_Tom20]

#'[foldchange > 2, log2foldchange > 1, qval 0.05  ]
dep <- add_rejections(data_diff_manual, alpha = 0.05, lfc = log2(1)) # but here it can only generate pval, padj

# Plot the first and second principal components
plot_pca(dep, x = 1, y = 2, n = 500, point_size = 4)

#View(plot_pca)

#'[Heatmap of all significant proteins]

# Plot a heatmap of all significant proteins with the data centered per protein
#'[with gene names] [each time same gene list]
a <- plot_heatmap(dep, type = "centered", kmeans = TRUE, 
                  k = 6, col_limit = 4, show_row_names = TRUE,
                  indicate = c("condition", "replicate"))

#'[without gene names]
plot_heatmap(dep, type = "centered", kmeans = TRUE, 
             k = 6, col_limit = 4, show_row_names = FALSE,
             indicate = c("condition", "replicate"))


#'[4.10.4Volcano plots of specific contrasts]

# Plot a volcano plot for the contrast 'Zuc_vs_Ctrl'
plot_volcano(dep, contrast = 'Zuc_vs_control', label_size = 4, add_names = FALSE)

plot_volcano(dep, contrast = 'Tom20_vs_control', label_size = 4, add_names = TRUE)

plot_volcano(dep, contrast = 'Zuc_vs_Tom20', label_size = 4, add_names = TRUE)

# Generate a results table
data_results <- get_results(dep)
# Generate a wide data.frame
df_wide <- get_df_wide(dep)

# write.csv(df_wide, "germ_line_cell_7062mata4_turboID_df_wide_deg_mass_spec_nes_zuc_tom20_MinDET_imputation_dep_20211125.csv")

#'[----------------------------------------------------------------]
#'[----------------------------------------------------------------]
getwd() #
library(EnhancedVolcano)
nes_vs_zuc_tom20_germcell_MinDet <- read.csv('germ_line_cell_7062mata4_turboID_df_wide_deg_mass_spec_nes_zuc_tom20_MinDET_imputation_dep_20211125.csv', header = TRUE)
colnames(nes_vs_zuc_tom20_germcell_MinDet)
# [1] "X.2"                          "name"                         "control_1"                    "control_2"                    "control_3"                    "Tom20_1"                     
# [7] "Tom20_2"                      "Tom20_3"                      "Zuc_1"                        "Zuc_2"                        "Zuc_3"                        "X.1"                         
# [13] "X"                            "Protein"                      "Proteinnames"                 "Entryname"                    "Genenames"                    "Entry.name"                  
# [19] "Status"                       "Protein.names"                "A"                            "B"                            "Organism"                     "Length"                      
# [25] "ID"                           "imputed"                      "num_NAs"                      "Tom20_vs_control_CI.L"        "Tom20_vs_control_CI.R"        "Tom20_vs_control_diff"       
# [31] "Tom20_vs_control_p.adj"       "Tom20_vs_control_p.val"       "Zuc_vs_control_CI.L"          "Zuc_vs_control_CI.R"          "Zuc_vs_control_diff"          "Zuc_vs_control_p.adj"        
# [37] "Zuc_vs_control_p.val"         "Zuc_vs_Tom20_CI.L"            "Zuc_vs_Tom20_CI.R"            "Zuc_vs_Tom20_diff"            "Zuc_vs_Tom20_p.adj"           "Zuc_vs_Tom20_p.val"          
# [43] "Tom20_vs_control_significant" "Zuc_vs_control_significant"   "Zuc_vs_Tom20_significant"     "significant"             
#''[-------------------------------------------------------------------]
library(qvalue)
#'#'[----------------------qvalue for nes_vs_zuc_germ_Cells_mindet---------------]
pvalues <- nes_vs_zuc_tom20_germcell_MinDet$Zuc_vs_control_p.val
qobj <- qvalue(p = pvalues)
nes_vs_zuc_tom20_germcell_MinDet$qval_Zuc_vs_control_germcell_minDet <- qobj$qvalues

library(tidyverse)
nes_vs_zuc_tom20_germcell_MinDet %>% select(, Zuc_vs_control_p.val, qval_Zuc_vs_control_germcell_minDet )

#     Zuc_vs_control_p.val qval_Zuc_vs_control_germcell_minDet
# 1           6.195285e-01                        1.760934e-01
# 2           3.333182e-02                        2.213063e-02
# 3           3.659866e-04                        7.234352e-04
#'[-----------------------------------------------------------]
#'[-----------------------------------------------------------]

#''[-------------------------------------------------------------------]
#'#'[----------------------qvalue for nes_vs_tom20_germ_Cells_mindet---------------]
pvalues <- nes_vs_zuc_tom20_germcell_MinDet$Tom20_vs_control_p.val
qobj <- qvalue(p = pvalues)
nes_vs_zuc_tom20_germcell_MinDet$qval_Tom20_vs_control_germcell_minDet <- qobj$qvalues

library(tidyverse)
nes_vs_zuc_tom20_germcell_MinDet %>% select(, Tom20_vs_control_p.val, qval_Tom20_vs_control_germcell_minDet )

#    Tom20_vs_control_p.val qval_Tom20_vs_control_germcell_minDet
#1             8.579614e-01                          2.209787e-01
#2             8.035571e-07                          1.819044e-06
#3             9.084874e-07                          1.917792e-06

#'[-----------------------------------------------------------]
#'[-----------------------------------------------------------]


#''[-------------------------------------------------------------------]
#'#'[----------------------qvalue for zuc_vs_tom20_germ_Cells_mindet---------------]
pvalues <- nes_vs_zuc_tom20_germcell_MinDet$Zuc_vs_Tom20_p.val
qobj <- qvalue(p = pvalues)
nes_vs_zuc_tom20_germcell_MinDet$qval_Zuc_vs_Tom20_germcell_minDet <- qobj$qvalues

library(tidyverse)
nes_vs_zuc_tom20_germcell_MinDet %>% select(, Zuc_vs_Tom20_p.val, qval_Zuc_vs_Tom20_germcell_minDet )
#    Zuc_vs_Tom20_p.val qval_Zuc_vs_Tom20_germcell_minDet
#1         7.486486e-01                      2.776596e-01
#2         3.586366e-06                      1.402667e-05
#3         5.443346e-05                      1.209405e-04
#     

# write.csv(nes_vs_zuc_tom20_germcell_MinDet, "germ_line_cell_7062mata4_turboID_df_wide_deg_mass_spec_nes_zuc_tom20_MinDET_imputation_dep_with_qvalue_20211125.csv")

#'[-----------------------------------------------------------]
#'[-----------------------------------------------------------]

























#'[-------------------------------------------------------------------]
#'[Figure 2]

library("dplyr")
library("DEP")
library("tidyr")
library("purrr")
library("ggplot2")
library("SummarizedExperiment")















#'[-------------------------------------------------------------------]
#'[Figure 2A]



#'[-------------------------------------------------------------------]
#'[Figure 2B]



#'[-------------------------------------------------------------------]
#'[Figure 2C]


#'[-------------------------------------------------------------------]
#'[-------------------------------------------------------------------]
#'[Figure 2D]

library(DEP)
library(ggplot2)
library(tidyverse)

nes_vs_zuc_tom20_germcell_MinDet <- read.csv("germ_line_cell_7062mata4_turboID_df_wide_deg_mass_spec_nes_zuc_tom20_MinDET_imputation_dep_with_qvalue_20211125.csv", header = TRUE)
colnames(nes_vs_zuc_tom20_germcell_MinDet)
# [1] "X.3"                                   "X.2"                                   "name"                                  "control_1"                             "control_2"                            
# [6] "control_3"                             "Tom20_1"                               "Tom20_2"                               "Tom20_3"                               "Zuc_1"                                
# [11] "Zuc_2"                                 "Zuc_3"                                 "X.1"                                   "X"                                     "Protein"                              
# [16] "Proteinnames"                          "Entryname"                             "Genenames"                             "Entry.name"                            "Status"                               
# [21] "Protein.names"                         "A"                                     "B"                                     "Organism"                              "Length"                               
# [26] "ID"                                    "imputed"                               "num_NAs"                               "Tom20_vs_control_CI.L"                 "Tom20_vs_control_CI.R"                
# [31] "Tom20_vs_control_diff"                 "Tom20_vs_control_p.adj"                "Tom20_vs_control_p.val"                "Zuc_vs_control_CI.L"                   "Zuc_vs_control_CI.R"                  
# [36] "Zuc_vs_control_diff"                   "Zuc_vs_control_p.adj"                  "Zuc_vs_control_p.val"                  "Zuc_vs_Tom20_CI.L"                     "Zuc_vs_Tom20_CI.R"                    
# [41] "Zuc_vs_Tom20_diff"                     "Zuc_vs_Tom20_p.adj"                    "Zuc_vs_Tom20_p.val"                    "Tom20_vs_control_significant"          "Zuc_vs_control_significant"           
# [46] "Zuc_vs_Tom20_significant"              "significant"                           "qval_Zuc_vs_control_germcell_minDet"   "qval_Tom20_vs_control_germcell_minDet" "qval_Zuc_vs_Tom20_germcell_minDet" 



g0 <- subset(nes_vs_zuc_tom20_germcell_MinDet, name == 'zuc')
g1 <- subset(nes_vs_zuc_tom20_germcell_MinDet, name == 'Gasz')
g2 <- subset(nes_vs_zuc_tom20_germcell_MinDet, name == 'mino')
g4 <- subset(nes_vs_zuc_tom20_germcell_MinDet, name == 'CG11513-RB') #'[Armi]
g5 <- subset(nes_vs_zuc_tom20_germcell_MinDet, name == 'SoYb')
g6 <- subset(nes_vs_zuc_tom20_germcell_MinDet, name == 'shu')
g16 <- subset(nes_vs_zuc_tom20_germcell_MinDet, name == 'daed')


library(DEP)
#'[https://stackoverflow.com/questions/58137828/create-new-variable-by-multiple-conditions-via-mutate-case-when]


p <- ggplot(nes_vs_zuc_tom20_germcell_MinDet, aes(Zuc_vs_control_diff, -log10(qval_Zuc_vs_control_germcell_minDet))) + geom_vline(xintercept = 0) + geom_point(colour = 'grey')+
  geom_point(aes(color = dplyr::case_when(qval_Zuc_vs_control_germcell_minDet < 0.05 & Zuc_vs_control_diff > 1 ~ "gray50", 
                                          qval_Zuc_vs_control_germcell_minDet < 0.05 & Zuc_vs_control_diff < -1 ~ "gray50",
                                          TRUE ~ "gray90")))  + theme_DEP1()  + ylim(0, 8) +  theme(legend.position = 'none') +
  scale_color_identity() + xlab('log2 fold-change') + ylab(bquote('-log10 qval')) +  geom_vline(xintercept= c(-1, 1), colour= 'black', linetype= 'dashed') + geom_hline(yintercept= -log10(0.05), colour= 'black', linetype= 'dashed') + ggtitle("Germ cell, MinDet, nes vs zuc, qval < 0.05, log2fc > 1") + annotate("text", x = 10, y=0, label = "Zuc") +  annotate("text", x = -7, y=0, label = "Nes") + 
  geom_point(data=g0, colour="red") +  # this adds a red point
  geom_text(data=g0, label="zuc", vjust=1, fontface = "bold") + # this adds a label for the red point
  geom_point(data=g1, colour="red") +  # this adds a red point
  geom_text(data=g1, label="Gasz", vjust=1, fontface = "bold") +
  
  geom_point(data=g2, colour="red") +  # this adds a red point
  geom_text(data=g2, label="mino", vjust=1, fontface = "bold") +
  
  geom_point(data=g4, colour="red") +  # this adds a red point
  geom_text(data=g4, label="Armi", fontface = "bold",vjust=1) +
  
  geom_point(data=g5, colour="red") +  # this adds a red point
  geom_text(data=g5, label="SoYb", vjust=1, fontface = "bold") +
  
  geom_point(data=g6, colour="red") +  # this adds a red point
  geom_text(data=g6, label="shu", vjust=1, fontface = "bold") +

  geom_point(data=g16, colour="red") +  # this adds a red point
  geom_text(data=g16, label="Daed", vjust=1, fontface = "bold") 

p

#'[-------------------------------------------------------------------]
#'[-------------------------------------------------------------------]
#'[Figure 2F]

#'[-------------------------------------------------------------------]
#'[Figure 2G]





#'[-------------------------------------------------------------------]
#'[-------------------------------------------------------------------]
#'[Figure 3A]

library(tidyverse)
library(ggplot2)
library(DEP)
nes_vs_zuc_tom20_germcell_MinDet <- read.csv("germ_line_cell_7062mata4_turboID_df_wide_deg_mass_spec_nes_zuc_tom20_MinDET_imputation_dep_with_qvalue_20211125.csv", header = TRUE)
colnames(nes_vs_zuc_tom20_germcell_MinDet)
# [1] "X.3"                                   "X.2"                                   "name"                                  "control_1"                             "control_2"                            
# [6] "control_3"                             "Tom20_1"                               "Tom20_2"                               "Tom20_3"                               "Zuc_1"                                
# [11] "Zuc_2"                                 "Zuc_3"                                 "X.1"                                   "X"                                     "Protein"                              
# [16] "Proteinnames"                          "Entryname"                             "Genenames"                             "Entry.name"                            "Status"                               
# [21] "Protein.names"                         "A"                                     "B"                                     "Organism"                              "Length"                               
# [26] "ID"                                    "imputed"                               "num_NAs"                               "Tom20_vs_control_CI.L"                 "Tom20_vs_control_CI.R"                
# [31] "Tom20_vs_control_diff"                 "Tom20_vs_control_p.adj"                "Tom20_vs_control_p.val"                "Zuc_vs_control_CI.L"                   "Zuc_vs_control_CI.R"                  
# [36] "Zuc_vs_control_diff"                   "Zuc_vs_control_p.adj"                  "Zuc_vs_control_p.val"                  "Zuc_vs_Tom20_CI.L"                     "Zuc_vs_Tom20_CI.R"                    
# [41] "Zuc_vs_Tom20_diff"                     "Zuc_vs_Tom20_p.adj"                    "Zuc_vs_Tom20_p.val"                    "Tom20_vs_control_significant"          "Zuc_vs_control_significant"           
# [46] "Zuc_vs_Tom20_significant"              "significant"                           "qval_Zuc_vs_control_germcell_minDet"   "qval_Tom20_vs_control_germcell_minDet" "qval_Zuc_vs_Tom20_germcell_minDet" 




g1 <- subset(nes_vs_zuc_tom20_germcell_MinDet, name == 'Marf')
g2 <- subset(nes_vs_zuc_tom20_germcell_MinDet, name == 'Tom20')
g3 <- subset(nes_vs_zuc_tom20_germcell_MinDet, name == 'spoon')
g4 <- subset(nes_vs_zuc_tom20_germcell_MinDet, name == 'porin')

g6 <- subset(nes_vs_zuc_tom20_germcell_MinDet, name == 'Mul1')
g7 <- subset(nes_vs_zuc_tom20_germcell_MinDet, name == 'Tom70')
g9 <- subset(nes_vs_zuc_tom20_germcell_MinDet, name == 'Tom40')

p <- ggplot(nes_vs_zuc_tom20_germcell_MinDet, aes(Tom20_vs_control_diff, -log10(qval_Tom20_vs_control_germcell_minDet))) + geom_vline(xintercept = 0) + geom_point(colour = 'grey')+
  geom_point(aes(color = dplyr::case_when(qval_Tom20_vs_control_germcell_minDet < 0.05 & Tom20_vs_control_diff > 1 ~ "gray50", 
                                          qval_Tom20_vs_control_germcell_minDet < 0.05 & Tom20_vs_control_diff < -1 ~ "gray50",
                                          TRUE ~ "gray90")))  + theme_DEP1()  + ylim(0, 8) + xlim(c(-8, 10))+ theme(legend.position = 'none') +
  scale_color_identity() + xlab('log2 fold-change') + ylab(bquote('-log10 qval')) +  geom_vline(xintercept= c(-1, 1), colour= 'black', linetype= 'dashed') + geom_hline(yintercept= -log10(0.05), colour= 'black', linetype= 'dashed') + ggtitle("Tom20 vs NES") + annotate("text", x = 10, y=0, label = "Tom20") +  annotate("text", x = -7, y=0, label = "Nes") + 

  geom_point(data=g1, colour="red") +  # this adds a red point
  geom_text(data=g1, label="Marf", vjust=1, fontface = "bold") +
  
  geom_point(data=g2, colour="red") +  # this adds a red point
  geom_text(data=g2, label="Tom20", vjust=1, fontface = "bold") +
  
  geom_point(data=g3, colour="red") +  # this adds a red point
  geom_text(data=g3, label="spoon", vjust=1, fontface = "bold") +
  
  geom_point(data=g4, colour="red") +  # this adds a red point
  geom_text(data=g4, label="Porin", fontface = "bold",vjust=1) +
  
  geom_point(data=g6, colour="red") +  # this adds a red point
  geom_text(data=g6, label="Mul1", vjust=1, fontface = "bold") +
  
  geom_point(data=g7, colour="red") +  # this adds a red point
  geom_text(data=g7, label="Tom70", vjust=1, fontface = "bold") +
  
  # geom_point(data=g8, colour="red") +  # this adds a red point
  # geom_text(data=g8, label="clu", vjust=1, fontface = "bold") +
  
  geom_point(data=g9, colour="red") +  # this adds a red point
  geom_text(data=g9, label="Tom40", vjust=1, fontface = "bold") 

p

#'[-------------------------------------------------------------------]
#'[-------------------------------------------------------------------]
#'[-------------------------------------------------------------------]
#'[-------------------------------------------------------------------]
#'[-------------------------------------------------------------------]
#'[-------------------------------------------------------------------]
#'[-------------------------------------------------------------------]
#'[-------------------------------------------------------------------]
#'[-------------------------------------------------------------------]
#'[-------------------------------------------------------------------]
#'[-------------------------------------------------------------------]
#'[-------------------------------------------------------------------]
#'[-------------------------------------------------------------------]
#'[-------------------------------------------------------------------]
#'[-------------------------------------------------------------------]
#'[-------------------------------------------------------------------]
#'[-------------------------------------------------------------------]
#'[-------------------------------------------------------------------]
#'[-------------------------------------------------------------------]
#'[-------------------------------------------------------------------]
#'[-------------------------------------------------------------------]











#'[-------------------------------------------------------------------]
#'[Figure 3C]

#'[-------------------------------------------------------------------]
#'[Figure 3D]


#'[-------------------------------------------------------------------]
#'[Figure 4A]

#'[-------------------------------------------------------------------]
#'[Figure 5C]


#'[-------------------------------------------------------------------]
#'[Figure 5F]


#'[-------------------------------------------------------------------]



#'[-------------------------------------------------------------------]
#'[Supplementary Figure 2B]

#'[-------------------------------------------------------------------]
#'[Supplementary Figure 2C]







#'[-------------------------------------------------------------------]
#'[-------------------------------------------------------------------]
#'[-------------------------------------------------------------------]























