#'  The following code pulls the MZmine data table, Sirius output, and metadata into usable dataframes
#'  to create graphs, correlations, etc.


# Call libraries needed for R
if (!require(dplyr)) install.packages('dplyr')
library(dplyr)
if (!require(tidyr)) install.packages('tidyr')
library(tidyr)
if (!require(readr)) install.packages('readr')
library(readr)
if (!require(stringr)) install.packages('stringr')
library(stringr)
if (!require(ecodist)) install.packages('ecodist')
library(ecodist)
if (!require(ggthemes)) install.packages('ggthemes')
library(ggthemes)
if (!require(ggplot2)) install.packages('ggplot2')
library(ggplot2)
if (!require(RColorBrewer)) install.packages('RColorBrewer')
library(RColorBrewer)
if (!require(svglite)) install.packages('svglite')
library(svglite)


####  1) Set working directory and pull files into R

setwd("~/Graduate School/Projects/CCE-LTER Cruise/P1706/Up to Date Code/P1706_index_ralph")

ASV16S <- read.csv('P1706_ASV16S.csv')
ASV18SV9 <- read.csv('P1706_ASV18SV9.csv')

ft <- read.csv('CCE1706_MZmine3_GNPS_1_quant.csv', header = T, check.names = F)
md <- read.csv('metadata_CCE.csv', header = T, check.names = F)

sirius <- read.csv('formula_identifications.tsv', sep = "\t", header = T, check.names = F)
canopus <- read.csv('canopus_formula_summary.tsv', sep = "\t", header = T, check.names = F)

#Arranging Sirius and Canopus files in the right format for further analysis:

#Renaming certain column headers to retain different molecular formulas
sirius <- rename(sirius, 'Sirius_molecularFormula' = 'molecularFormula')
canopus <- rename(canopus, 'Canopus_molecularFormula' = 'molecularFormula')

#Merging two Canopus and Sirius dataframes by id
Siriusall <- merge(canopus, sirius, by = 'id', all = TRUE)

# Indexing featureID out of id
Siriusall <- separate(data = Siriusall, col = id, into = c(NA, NA, NA, NA, NA, 'row.ID'), sep = '_')

# Indexing molecular formulas into C, N, H, P, O, S, etc.

# First, Sirius molecular formulas
molecules <- regmatches(Siriusall$Sirius_molecularFormula, gregexpr("\\b[A-Z][a-z]*\\d*", Siriusall$Sirius_molecularFormula))
molecules <- lapply(molecules, function(a) paste0(a, ifelse(grepl("[^0-9]$", a), "1", "")))

atomcounts <- lapply(molecules, function(mol) setNames(as.integer(gsub("\\D", "", mol)), gsub("\\d", "", mol)))

atoms <- unique(unlist(sapply(atomcounts, names)))
atoms <- sapply(atoms, function(atom) sapply(atomcounts, function(a) if (atom %in% names(a)) a[atom] else 0))
rownames(atoms) <- Siriusall$Sirius_molecularFormula
atoms

Sirius_formula <- as.data.frame(atoms)
rownames(Sirius_formula) <- NULL
colnames(Sirius_formula) <- paste(colnames(Sirius_formula), "Sirius", sep = '_')

Sirius_final <- cbind(Siriusall, Sirius_formula)

#Next, Canopus molecular formulas
molecules <- regmatches(Sirius_final$Canopus_molecularFormula, gregexpr("\\b[A-Z][a-z]*\\d*", Sirius_final$Canopus_molecularFormula))
molecules <- lapply(molecules, function(a) paste0(a, ifelse(grepl("[^0-9]$", a), "1", "")))

atomcounts <- lapply(molecules, function(mol) setNames(as.integer(gsub("\\D", "", mol)), gsub("\\d", "", mol)))

atoms <- unique(unlist(sapply(atomcounts, names)))
atoms <- sapply(atoms, function(atom) sapply(atomcounts, function(a) if (atom %in% names(a)) a[atom] else 0))
rownames(atoms) <- Sirius_final$Canopus_molecularFormula
atoms

Canopus_formula <- as.data.frame(atoms)
rownames(Canopus_formula) <- NULL
colnames(Canopus_formula) <- paste(colnames(Canopus_formula), "Canopus", sep = '_')

Sirius_final <- cbind(Sirius_final, Canopus_formula)

#Sirius_final is best I have currently

#________________________________________________________________________________________________________________________________________________________

#Van Krevelen Plots

df <- Sirius_final


#C:N, O:C, H:C and average C oxidation columns
df$C_N <- df$C_Sirius / df$N_Sirius
df$O_C <- df$O_Sirius / df$C_Sirius
df$H_C <- df$H_Sirius / df$C_Sirius
df$avCox <- -(1*df$H_Sirius - 3*df$N_Sirius - 2*df$O_Sirius + 5*df$P_Sirius -2*df$S_Sirius)/df$C_Sirius

#Only using Sirius formulas for cutoff, creating van krevelen data frame and cutoffs.
#Cutoff will be C/N < 4 ; O/C > 2; H/C > 2.5
df_vk <- subset(df, O_C <2 & H_C < 2.5 & C_N > 4.5)
df_vk <- subset(df_vk, avCox >-4)
df_vk <- subset(df_vk, avCox <4)

## Van Krevelen Plots
vk <- df_vk[-c(307:353, 356)]

#Molecular Formula
VK <- ggplot(vk, aes(x = O_C, y = H_C)) +
  geom_point(aes(color = ionMass), size = 2, na.rm = TRUE, alpha = 0.6) +
  scale_color_gradient(name = "Molecular Formula",
                       low = ("orange"), high = ("black")) +
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.title = element_text(face = "bold"), 
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), 
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24)) +
  ggtitle('Van Krevelen Plot of Sirius Formulas (4,931)') +
  labs(x = "O/C", y = "H/C") +
  scale_x_continuous(limits = c(0, 1.5), breaks = seq(0.0, 1.2, by = 0.3)) +
  scale_y_continuous(limits = c(0, 2.5), breaks = seq(0.0, 2.5, by = 0.5))
VK

#Avg Oxidation State

VK2 <- ggplot(df_vk, aes(x = O_C, y = H_C)) +
  geom_point(aes(color = avCox), size = 2, na.rm = TRUE, alpha = 0.6) +
  scale_color_gradient2(name = "Carbon Oxidation State",
                        low = ("blue"), mid = ("white"), high = ("red"), midpoint = 0) +
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.title = element_text(face = "bold"), 
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), 
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24)) +
  ggtitle('Van Krevelen Plot of Sirius Formulas (4,931)') +
  labs(x = "O/C", y = "H/C") +
  scale_x_continuous(limits = c(0, 1.5), breaks = seq(0.0, 1.2, by = 0.3)) +
  scale_y_continuous(limits = c(0, 2.5), breaks = seq(0.0, 2.5, by = 0.5))
VK2

#________________________________________________________________________________________________________________________________________________________

##################################ABZER CODE START#################################
# Arranging the metadata table in the right format for further analysis:

md <- md[, colSums(is.na(md))<nrow(md)] #removing NA columns in md, if any present
rownames(md) <- md$filename
md <- md[, -1]
rownames(md) <- gsub('Blank_', '', rownames(md)) #Removing "Blank" in the 1st two row names of md (Specific to dataset)


# Arranging the feature table in the right format for further analysis (Specific to dataset)

new_ft <- ft[,grep('mzXML',colnames(ft))] # only picking mzXML files
new_ft <- new_ft[,-grep('_STKL',colnames(new_ft))] # excluding the STKL columns
colnames(new_ft) <- gsub('_MSMS.mzXML Peak area','.mzxml',colnames(new_ft))  #substituting the file extension in the colnames of new_ft with mzxml
rownames(new_ft) <- paste(ft$'row ID',round(ft$'row m/z',digits = 3),round(ft$'row retention time',digits = 3), sep = '_') #Taking row ID, m/z value and RT of ft as the rownames of new_ft

new_ft <- new_ft[,order(colnames(new_ft))] #ordering the columns by names
new_ft <- new_ft[,-1] # Removing the column "Brandon_Deep_DOM"

#Merging feature table with metadata file

ft_final <- new_ft[,which(colnames(new_ft)%in%rownames(md))] #picking only the columns in new_ft that are comparable to rownames of md
ft_final  <- new_ft[,match(rownames(md),colnames(new_ft))] #returns the matched position of rownames of md to that of colnames of new_ft
identical(colnames(ft_final),rownames(md)) #checking if the colnames of ft_final and rownames of md matched. Should return TRUE
dim(ft_final) #column number should match below row number
dim(md) #row number should match above column number

#Creating frequencyPlot function to show distribution of frequency

#'Global' settings for plot size in the output cell
options(repr.plot.width=5, repr.plot.height=3)
FrequencyPlot <- function(x1,x2){
  #creating bins from -1 to 10^10 using sequence function seq()
  bins <- c(-1,0,(1 * 10^(seq(0,10,1)))) 
  #cut function cuts the give table into its appropriate bins
  scores_x1 <- cut(as.matrix(x1),bins,labels = c('0','1','10','1E2','1E3','1E4','1E5','1E6','1E7','1E8','1E9','1E10')) 
  #transform function convert the tables into a column format: easy for visualization 
  Table_x1<-transform(table(scores_x1)) #contains 2 columns: "scores_x1", "Freq"
  #Repeating the same steps for x2
  scores_x2 <- cut(as.matrix(x2),bins,labels = c('0','1','10','1E2','1E3','1E4','1E5','1E6','1E7','1E8','1E9','1E10'))
  Table_x2<-transform(table(scores_x2))
  #Getting the names of x1 and x2
  arg1 <- deparse(substitute(x1))
  arg2 <-deparse(substitute(x2))
  #Creating a data frame for plotting
  data_plot <- as.data.frame(c(Table_x1$Freq,Table_x2$Freq)) #Concatenating the frequency info of both tables rowwise
  colnames(data_plot) <- "Freq" #naming the 1st column as 'Freq'
  data_plot$Condition <- c(rep(arg1,12),rep(arg2,12)) #adding a 2nd column 'Condition', which just repeats the name of x1 and x2 accordingly
  data_plot$Range_bins <- rep(Table_x1$scores_x1,2) #Adding 3rd column 'Range Bins'
  data_plot$Log_Freq <- log(data_plot$Freq+1) #Log scaling the frequency values
  ## GGPLOT2
  BarPlot <- ggplot(data_plot, aes(Range_bins, Log_Freq, fill = Condition)) + 
    geom_bar(stat="identity", position = "dodge", width=0.4) + 
    scale_fill_brewer(palette = "Set1") +
    ggtitle(label="Frequency plot") +
    xlab("Range") + ylab("(Log)Frequency") + labs(fill = "Data Type") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +   # setting the angle for the x label
    theme(axis.text.y = element_text(angle = 45, vjust = 0.5, hjust=1)) +   # setting the angle for the y label
    theme(plot.title = element_text(hjust = 0.5)) # centering the plot title
  print(BarPlot)
}

#________________________________________________________________________________________________________________________________________________________

#' Blank Removal using Blank (Control) Samples
# A common filtering method is to use a cutoff to remove features that are not present sufficient enough in biological samples.
# 
# 1) We find an average for all the feature intensities in your control set and sample set. Therefore, for n no.of features in 
#    a control or sample set, we get n no.of averaged features.
# 2) Next, we get a ratio of this average_control vs average_sample. This ratio Control/sample tells us how much of that particular 
#    feature of a sample gets its contribution from control. If it is more than 30% (or Cutoff as 0.3), we consider the feature as 
#    noise/contamination.
# 3) The resultant information (if ratio > Cutoff or not) is stored in a bin
# 4) We count the no.of features in the bin that satisfies the condition ratio > cutoff, and consider those features as 'noise or background features'.
# 5) We also try to visualize the frequency distribution in our data Ctrl and samples

#Splitting the samples: 
Ctrl <-  ft_final[,1:2]
Samples <- ft_final[,-1:-2]

#Select Y to perform Blank Correction, Cutoff - 0.3 represents 30% cutoff, 
input_data <- ft_final

if(readline('Do you want to perform Blank Removal Process - Y/N:')=='Y'){
  Cutoff <- as.numeric(readline('Enter Cutoff value between 0.1 & 1:')) # (i.e. 10% - 100%). Ideal cutoff range: 0.1-0.3
  #When cutoff is low, more noise (or background) detected; With higher cutoff, less background detected, thus more features observed
  
  Avg_ctrl <- rowMeans(Ctrl, na.rm= FALSE, dims = 1) # set na.rm = FALSE to check if there are NA values. Because when set as TRUE, NA values are changed to 0
  Avg_samples <- rowMeans(Samples, na.rm= FALSE, dims = 1)
  Ratio_Ctrl_Sample <- (Avg_ctrl+1)/(Avg_samples+1)
  Bg_bin <- ifelse(Ratio_Ctrl_Sample > Cutoff, 1, 0 )  # checks if the Ratio is greater than Cutoff, if so put 1, else 0 in Bg_bin
  
  # to check if there are any NA values present. It is not good to have NA values in the 4 variables as it will affect the final dataset to be created
  temp_NA_Count <-cbind(Avg_ctrl ,Avg_samples,Ratio_Ctrl_Sample,Bg_bin)
  
  print('No of NA values in the following columns:')
  print(colSums(is.na(temp_NA_Count)))
  
  Bg_Features <- sum(Bg_bin ==1,na.rm = TRUE) # Calculates the number of background features present
  No_of_Features <- nrow(input_data) - Bg_Features
  
  print(paste("No.of Background or noise features:",Bg_Features))
  print(paste("No.of features after excluding noise:",No_of_Features)) 
  
  input_data1 <- cbind(as.matrix(input_data),Bg_bin)    
  plot_CtrlSample <- FrequencyPlot(Samples,Ctrl)
}



#This plot will show you where the LOD is (what is the minimum (other than 0) signal value shown?)
#________________________________________________________________________________________________________________________________________________________

#Imputation:
# For several reasons, real world datasets might have some missing values in it, in the form of NA, NANs or 0s. 
# Even though the gapfillig step of MZmine fills the missing values, we still end up with some missing values or 
# Os in our feature table. This could be problematic for statistical analysis. In order to have a better dataset, 
# we cannot simply discard those rows or columns with missing values as we will lose a chunk of our valuable data. 
# Instead we can try imputing those missing values. Imputation involves replacing the missing values in the data 
# with a meaningful, reasonable guess. There are several methods, such as:
#  Mean imputation (replacing the missing values in a column with the mean or average of the column)
#  Replacing it with the most frequent value
#  Several other machine learning imputation methods such as k-nearest neighbors algorithm(k-NN), Hidden Markov Model(HMM)
# The method that we use is: Replacing the zeros from the gapfilled quant table with the Cutoff_LOD value.

#Continuing from previous steps in Blank Removal :
#  1) From the plot, we decide on a cutoff_LOD value (LOD-Limit of Detection). If until range 10-100, (shown in the 
#     figure as 1E2) there are no or very less features, we want to exclude until that range and consider from range 
#     (100-1000), or, in other words, we take 'range:100-1000 or 1E3 or 1000' as Cutoff_LOD. This value will be used 
#     to replace the zeros in the data table
#  2) Once we consider that LOD value,
#      a) We create a temparory dataset with all the feature intensites of your sample (not the ctrl, only the sample) 
#      and checking it against the cutoff_LOD value. If it is less than the cutoff_LOD, we replace it with cutoff_LOD. 
#      Thus, for instance, if we take 1000 as cutoff_LOD, our sample data will be filled with a bunch of 1000s now instead of zeros
#      b) Then we create a Final dataset using the temporary dataset. Here, we try to see if each feature from all samples is noise or
#     not (from step 4), if it noise, we replace the feature with cutoff_LOD, else we keep the info from the temporary dataset.
# For the code below, I selected 'N' and gave a LOD value: 1000

#Enter the LOD value as seen in the frequency plot: Ex: 1E3 or 1000
Cutoff_LOD <-ifelse(readline("Was Imputation step already performed? Y/N")=="Y",RawLOD,as.numeric(readline("Enter your Cutoff LOD here:")))  
#________________________________________________________________________________________________________________________________________________________
# Replacing the Sample intensities with Cutoff_LOD, if they are lower than Cutoff_LOD
temp_matrix <- c()
for (i in 1:ncol(Samples)){ 
  x <- ifelse(Samples[,i] > Cutoff_LOD, Samples[,i],Cutoff_LOD)
  temp_matrix <- cbind(temp_matrix,as.matrix(x))
}
colnames(temp_matrix) <- colnames(Samples)

# Replacing the Sample intensities with Cutoff_LOD, if they are considered as a noise(bg_bin==1)
Final_matrix <-c()
for (i in 1:ncol(temp_matrix)){
  x <-ifelse(input_data1[,ncol(input_data1)] ==1, Cutoff_LOD, temp_matrix[,i])
  Final_matrix <-cbind(Final_matrix,x)
}
colnames(Final_matrix) <- colnames(Samples)

# Final_Matrix contains datafile that has undergone imputation

#Normalization:
  #The following code performs sample-centric (column-wise) normalisation:
  
  #Sample centric normalisation:----------------------------------------------- I said "Y" (YES) here
  if (readline("Do you want to perform Normalization: Y/N:") == 'Y'){  
    sample_sum <- colSums(Final_matrix, na.rm= TRUE, dims = 1)
    Normalized_data <- c()
    for (i in 1:ncol(Final_matrix)){
      x <- Final_matrix[,i] / sample_sum[i]
      Normalized_data <- cbind(Normalized_data, x)
    }
    
    colnames(Normalized_data) <- names(sample_sum)
    
    
  } else return(Final_matrix)
#________________________________________________________________________________________________________________________________________________________
print(paste('No.of NA values in Normalized data:',sum(is.na(Normalized_data)== TRUE)))
#write.csv(Normalized_data,file=paste0('CCE_Normalised_Data' ,'.csv'),row.names =FALSE)


## Multi-Dimensional Scaling (MDS) Plots
#Visualize multivariate data in 2D plots

#excluding Blanks from md and some other filtering conditions
MetaData <-  md %>% filter(ATTRIBUTE_Filament_Possition != 'Blank',
                           ATTRIBUTE_Location!='Santa_Barbara_Basin',
                           ATTRIBUTE_Location!='Transcet_1',
                           ATTRIBUTE_Location!='Transcet_2',
                           ATTRIBUTE_Location!='Transect_3',
                           ATTRIBUTE_Depth <= 20 ) 

# the corresponding column files for the filtered MetaData is picked from the normalized data
md_data <- Normalized_data[,which(colnames(Normalized_data)%in%rownames(MetaData))]

#For PCA plots, use Euclidean distance:
#dist_matrix <- dist(t(md_data), method="euclidean")

#For PCoA plots using Bray-curtis distance:
dist_matrix <- bcdist(t(md_data)) # transposed in order to compute the distance between the columns of a data matrix

pcoa<- cmdscale(dist_matrix, eig = TRUE, x.ret=TRUE)
pcoa.var.per <-round(pcoa$eig/sum(pcoa$eig)*100,1)
pcoa.values <- pcoa$points
pcoa.data <- data.frame(MetaData$ATTRIBUTE_Depth_Range,
                        X=pcoa.values[,1],
                        Y=pcoa.values[,2])

##GGPLOT
Plot <- ggplot(pcoa.data, aes(x=X, y=Y, col= as.factor(MetaData$ATTRIBUTE_Filament_Possition))) + 
  #geom_jitter(aes(shape = as.factor(MetaData$ATTRIBUTE_Depth)), size = 3) +
  geom_point(size=2,alpha=0.8)  +
  ggtitle(label="MDS plot") +
  xlab(paste0("MDS1 : ",pcoa.var.per[1],"%",sep="")) + 
  ylab(paste0("MDS2 : ",pcoa.var.per[2],"%",sep="")) + 
  labs(color = "Cycles",shape='Depth Range') + 
  theme(plot.title = element_text(hjust = 0.5)) 

Plot <- Plot + labs(subtitle = 'Using Bray_Curtis Distance on normalized data')
Plot
#ggsave(file="PCoA_Depth_0_20.svg", plot=Plot, width=10, height=8)

#To visualize the color palette

#options(repr.plot.width=5, repr.plot.height=6) #'global' settings for plot size in the output cell
#display.brewer.all()

# In the below code, n=5 is selected as default. Because, in general, the palette breaks the color 
# from darker shade to lighter. In order to avoid considering really light colors in the plot (as 
# it is hard to see), we specify a random n=5 as it is greater than no.of days (Max no.of days in 
# a cycle=4) and take the first few colors according to the no.of days in each cycle

a1 <- rev(brewer.pal(5, "YlOrRd"))
b1 <- rev(brewer.pal(5, "PuBu"))
c1 <- rev(brewer.pal(5, "Greens"))
d1 <- rev(brewer.pal(5, "RdPu"))

plot(rep(1,length(b1)),col= b1,pch=19,cex=3,ylab=NULL,xlab=NULL) #For visualising the gradient colors in a plot

Primary_level <- levels(factor(MetaData$ATTRIBUTE_Filament_Possition)) # Getting the levels "Cycle_1" to "Cycle_4"
print(Primary_level)
col_used <-  factor(MetaData$ATTRIBUTE_Location) # has the information of different days in each cycle

a1 <- a1[1:length(grep(Primary_level[1],levels(col_used)))]
b1 <- b1[1:length(grep(Primary_level[2],levels(col_used)))]
c1 <- c1[1:length(grep(Primary_level[3],levels(col_used)))]
d1 <- d1[1:length(grep(Primary_level[4],levels(col_used)))]

colors_gradient <- c(a1,b1,c1,d1)
colors_gradient

plot(rep(1,length(b1)),col= b1,pch=19,cex=3,ylab=NULL,xlab=NULL) #For visualising the gradient colors in a plot

Plot2 <- ggplot(pcoa.data, aes(x=X, y=Y, col= as.factor(MetaData$ATTRIBUTE_Location))) + 
  geom_point(size=2,alpha=0.8)  +
  ggtitle(label="MDS plot") +
  xlab(paste0("MDS1 : ",pcoa.var.per[1],"%",sep="")) + 
  ylab(paste0("MDS2 : ",pcoa.var.per[2],"%",sep="")) + 
  labs(color = "Location") + 
  theme(plot.title = element_text(hjust = 0.5)) 

Plot2 <- Plot2 + scale_color_manual(values = colors_gradient)+
  labs(subtitle = 'Using Bray_Curtis Distance on normalized data')

Plot2

#To save a ggplot as svg (scaled vector graph)
ggsave(file="PCoA_Depth_0_20.svg", plot=Plot2, width=10, height=8)

