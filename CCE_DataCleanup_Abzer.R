library(dplyr)
library(ggplot2)
library(ecodist) #For PcoA using Bray-Curtis distance


setwd(normalizePath('D:\Projects\CCE_Data-Analysis-main\P1706_index_ralph'),"/",mustWork=FALSE)
setwd('D:/Projects/CCE_Data-Analysis-main/P1706_index_ralph')


ft <- read.csv("CCE1706_MZmine3_GNPS_1_quant.csv",check.names = FALSE)
md <- read.csv('metadata_CCE.csv', check.names = FALSE,header = TRUE)


md <- md[,colSums(is.na(md))<nrow(md)] #removing NA columns in md, if any present
rownames(md) <- md$filename
md <- md[,-1]
rownames(md) <- gsub('Blank_','',rownames(md))

new_ft <- ft[,grep('mzXML',colnames(ft))] # only picking mzXML columns
new_ft <- new_ft[,-grep('_STKL',colnames(new_ft))] # excluding the STKL columns
colnames(new_ft) <- gsub('mzXML','mzxml',colnames(new_ft))
colnames(new_ft) <- gsub('_MSMS','',colnames(new_ft)) #removing _MSMS from the colnames of new_ft
rownames(new_ft) <- paste(ft$'row ID',round(ft$'row m/z',digits = 3),round(ft$'row retention time',digits = 3), sep = '_')

new_ft <- new_ft[,order(colnames(new_ft))] #ordering the columns by names
new_ft <- new_ft[,-1] # Removing the column "Brandon_Deep_DOM"

new_md <- md[which(rownames(md)%in%colnames(new_ft)),] #picking only the metadata rows by its rownames that are comparable to colnames of new_ft
ft_final <- new_ft[,which(colnames(new_ft)%in%rownames(md))] #picking only the columns in new_ft that are comparable to rownames of md
ft_final  <- new_ft[,match(rownames(md),colnames(new_ft))]
identical(colnames(ft_final),rownames(md)) #checking if the colnames of ft_final and rownames of md matched. Should return TRUE


#Splitting the samples----------
Ctrl <-  ft_final[,1:2]
Samples <- ft_final[,-1:-2]

##------------------ Frequency Plot Function------------------------

#The below function takes in the two input datatables.
#calculates the frequency distribution of the data in the order of 10
#and produces a grouped barplot showing the distribution as output. 
#The frequency plot shows where the features are present in higher number.

options(repr.plot.width=5, repr.plot.height=3) #'global' settings for plot size in the output cell

FrequencyPlot <- function(x1,x2){
  
  bins <- c(-1,0,(1 * 10^(seq(0,10,1)))) #creating bins from -1 to 10^10 using sequence function: seq()
  
  scores_x1 <- cut(as.matrix(x1),bins,labels = c('0','1','10','1E2','1E3','1E4','1E5','1E6','1E7','1E8','1E9','1E10')) #cut function store the data into the appropriate bins
  Table_x1<-transform(table(scores_x1)) #transform function convert the tables into column format: easy for visualization
  
  # Repeating the same steps for x2
  scores_x2 <- cut(as.matrix(x2),bins,labels = c('0','1','10','1E2','1E3','1E4','1E5','1E6','1E7','1E8','1E9','1E10'))
  Table_x2<-transform(table(scores_x2))
  
  arg1 <- deparse(substitute(x1))
  arg2 <-deparse(substitute(x2))
  
  data_plot <- as.data.frame(c(Table_x1$Freq,Table_x2$Freq))
  colnames(data_plot) <- "Freq"
  data_plot$Condition <- c(rep(arg1,12),rep(arg2,12))
  data_plot$Range_bins <- rep(Table_x1$scores_x1,2)
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


#--Imputation---------------------------------------------------------------------
GapFilled <- ft_final

if(readline('Imputation - YES or NO:')=='YES'){
  
  plot<- FrequencyPlot(GapFilled,NotGapFilled)
  
  Arg1 = plot$data$Condition[1]
  Arg2 = plot$data$Condition[13]
  plotData_New <- subset(plot$data,plot$data$Freq!=0 & plot$data$Range_bins !=0) # accessing the datatable of plot and subsetting with the condition: Eliminating the Range (or bin) 0 and Ranges with zero frequencies 
  First_val_temp <- aggregate(plotData_New$Freq, by=list(plotData_New$Condition), FUN=first) #getting the first appearing value of this new plot datatable
  First_val <- plotData_New[plotData_New$Freq %in% c(First_val_temp$x[1],First_val_temp$x[2]),] # Subsetting the rows in the plotData_New that has the first appearing values
  
  print(paste0("The Range with a minimum value greater than 0 for ",Arg1,":", First_val$Range_bins[1]))
  print(paste0("The Range with a minimum value greater than 0 for ",Arg2,":", First_val$Range_bins[2]))
  
  RawLOD <- round(min(NotGapFilled[NotGapFilled!=min(NotGapFilled)])) # getting the 2nd minimum value of non-gap filled data. (The first minimum value in the data table is usually zero)
  
  print(paste0("The minimum value in the Non-gap filled data other than 0 : ",RawLOD))
  GapFilled[GapFilled==0 | GapFilled<RawLOD] <- RawLOD # Replacing zeros in the gap-filled file as well as values<RawLOD with RawLOD
  RawLOD_Table <- GapFilled
  #write.csv(RawLOD_Table, file=paste0('Quant_Table_filled_with_MinValue_',RawLOD,'_NotGapFilled','.csv'),row.names =FALSE) 
  input_data <- GapFilled
} else input_data <- GapFilled




#---Blank Removal-------------------------------------------------------------------

if(readline('Blank Removal Process - YES or NO:')=='YES'){
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

#The frequency plot shows where the features are present in higher number. 
#For ex: if until 1E2 has no or really less features, the goal is to exclude until that range and consider only values from 1E3 range. 
#Thus 1E3 will be taken as Cutoff_LOD (Limit of Detection). This value will eventually become the 'new zero'.

#The following step asks if the imputation was already performed, if so, it takes that value as the Cutoff_LOD, 
#else, we get to specify our Cutoff_LOD based on the frequency plot. 
#Once we have our Cutoff_LOD: 
#--- We create a temporary dataset checking all the feature intensites of our sample (only the sample, without control) and check it against the Cutoff_LOD. 
#    If it is less than the Cutoff_LOD, we replace it with Cutoff_LOD. Thus our sample data, for example, is with a bunch of 1000s (if our Cutoff_LOD=1000) instead of zeros 
#--- Then, we create a Final dataset using the temporary dataset. 
#   Here, we try to see if each feature from all samples is noise or not. 
#   If it's noise, we replace it with Cutoff_LOD as well, else we keep the info from the temporary dataset as such.

{
  Cutoff_LOD <-ifelse(readline("Was Imputation step already performed?")=="YES",RawLOD,as.numeric(readline("Enter your Cutoff LOD here:")))  #Enter the LOD value as seen in the frequency plot
  temp_matrix <- c()
  for (i in 1:ncol(Samples)){ # Replacing the Sample intensities with Cutoff_LOD, if they are lower than Cutoff_LOD
    x <- ifelse(Samples[,i] > Cutoff_LOD, Samples[,i],Cutoff_LOD)
    temp_matrix <- cbind(temp_matrix,as.matrix(x))
  }
  colnames(temp_matrix) <- colnames(Samples)
  
  Final_matrix <-c()
  for (i in 1:ncol(temp_matrix)){
    x <-ifelse(input_data1[,ncol(input_data1)] ==1, Cutoff_LOD, temp_matrix[,i])
    Final_matrix <-cbind(Final_matrix,x)
  }
  colnames(Final_matrix) <- colnames(Samples)
}
write.csv(Final_matrix,file=paste0('Processed_Quant_Table_filled_with_Value_',Cutoff_LOD ,'.csv'),row.names =FALSE)



#Sample centric normalisation:-----------------------------------------------
if (readline("Normalization: YES or NO:") == 'YES'){  
  sample_sum <- colSums(Final_matrix, na.rm= TRUE, dims = 1)
  Normalized_data <- c()
  for (i in 1:ncol(Final_matrix)){
    x <- Final_matrix[,i] / sample_sum[i]
    Normalized_data <- cbind(Normalized_data, x)
  }
  
  colnames(Normalized_data) <- names(sample_sum)
  
  
} else return(Final_matrix)

print(paste('No.of NA values in Normalized data:',sum(is.na(Normalized_data)== TRUE)))


##--------PCoA plot: Using Bray-Curtis Distance------------------------------------------------------------------------------------------


MetaData <-  md %>% filter(ATTRIBUTE_Filament_Possition != 'Blank',ATTRIBUTE_Location!='Santa_Barbara_Basin',ATTRIBUTE_Location!='Transcet_1',ATTRIBUTE_Location!='Transcet_2',ATTRIBUTE_Location!='Transect_3') #excluding Blanks from md and some other filtering conditions
md_data <- Final_matrix[,which(colnames(Final_matrix)%in%rownames(MetaData))] # the corresponding column files for the filtered metadata is picked from the normalized data

dist_matrix <- bcdist(t(md_data)) # transposed in order to compute the distance between the columns of a data matrix
pcoa<- cmdscale(dist_matrix, eig = TRUE, x.ret=TRUE)
pcoa.var.per <-round(pcoa$eig/sum(pcoa$eig)*100,1)
pcoa.values <- pcoa$points
pcoa.data <- data.frame(MetaData$ATTRIBUTE_Location,
                        X=pcoa.values[,1],
                        Y=pcoa.values[,2])


Plot <- ggplot(pcoa.data, aes(x=X, y=Y, col= as.factor(MetaData$ATTRIBUTE_Location))) + 
  geom_point(size=4,alpha=0.8)  +
  ggtitle(label="MDS plot using Bray-Curtis Distance") +
  xlab(paste0("MDS1 : ",pcoa.var.per[1],"%",sep="")) + 
  ylab(paste0("MDS2 : ",pcoa.var.per[2],"%",sep="")) + 
  labs(color = "Location") + 
  theme(plot.title = element_text(hjust = 0.5)) 

Plot
