#Function: InsideLevels
 InsideLevels <- function(metatable){
    lev <- c()
    typ<-c()
    for(i in 1:ncol(metatable)){
      x <- levels(droplevels(as.factor(metatable[,i])))
      if(is.double(metatable[,i])==T){x=round(as.double(x),2)}
      x <-toString(x)
      lev <- rbind(lev,x)
      
      y <- class(metatable[,i])
      typ <- rbind(typ,y)
    }
    out <- data.frame(INDEX=c(1:ncol(metatable)),ATTRIBUTES=colnames(metatable),LEVELS=lev,TYPE=typ,row.names=NULL)
    return(out)
  }

#Function: SubsetLevels  
SubsetLevels <- function(metatable){
    IRdisplay::display(InsideLevels(metatable)) # use show() when working in RStudio
    Condition <- as.double(unlist(strsplit(readline("Enter the IDs of interested attributes to subset (separaed by commas if more than one attribute):"), split=',')))
    
    for( i in 1:length(Condition)){
      list_final <- c() 
      #Shows the different levels within each selected condition:
      Levels_Cdtn <- levels(droplevels(as.factor(metatable[,Condition[i]])))
      subset_meta <- data.frame(Index=1:length(Levels_Cdtn),Levels_Cdtn)
      colnames(subset_meta)[2] <- paste("Levels_",colnames(metatable[Condition[i]]))
      IRdisplay::display(subset_meta) # use show() when working in RStudio
      
      
      #Among the shown levels of an attribute, select the ones to keep or exclude:
      Read_cdtn <- readline("Do you want to keep or exclude few conditions? K/E: ")
      if( Read_cdtn=="K"){
        temp <- as.double(unlist(strsplit(readline("Enter the index numbers of condition(s) you want to KEEP (separated by commas):"), split=',')))
        ty <- class(metatable[,Condition[i]])
        list_keep <-Levels_Cdtn[temp, drop=F]
        list_exc <-Levels_Cdtn[-temp, drop=F]
        cat("The condition(s) you want to exclude in ",colnames(metatable)[Condition[i]]," :",list_exc,"\n")
        cat("The condition(s) you want to keep in ",colnames(metatable)[Condition[i]]," :",list_keep,"\n")
        
        }else if(Read_cdtn=="E"){
          temp <- as.double(unlist(strsplit(readline("Enter the index numbers of condition(s) you want to EXCLUDE (separated by commas):"), split=',')))
          ty <- class(metatable[,Condition[i]])
          list_exc <-Levels_Cdtn[temp, drop=F]
          list_keep <-Levels_Cdtn[-temp, drop=F]
          cat("The condition(s) you want to exclude in ",colnames(metatable)[Condition[i]]," :",list_exc,"\n")
          cat("The condition(s) you want to keep in ",colnames(metatable)[Condition[i]]," :",list_keep,"\n")
        
          }else{
            print("Sorry! You have given a wrong input!! Please enter either K or E")
            break
      }
      
      #In order to keep the original datatype of the columns, we define the following conditions, else it would all become 'characters or factors'
      if(ty=="integer"){list_keep <- as.integer(list_keep);list_exc <- as.integer(list_exc)
      }else if(ty=="double"){list_keep <- as.double(list_keep);list_exc <- as.double(list_exc)
      }else if(ty=="numeric"){list_keep <- as.numeric(list_keep);list_exc <- as.numeric(list_exc)}
      
      #Gets all the elements in list_keep into list_final 
      for(j in 1:length(list_keep)){
        sub_list <- metatable[(metatable[,Condition[i]] == list_keep[j]),]
        list_final <- rbind(list_final,sub_list)
      }
      metatable <- list_final #list_final again called as metatable in order to keep it in the for-loop for further subsetting
    }
    return(metatable)
}


#Function: FrequencyPlot

FrequencyPlot <- function(x1,x2){
  
  #creating bins from -1 to 10^10 using sequence function seq()
  bins <- c(-1,0,(1 * 10^(seq(0,10,1)))) 
  
  #cut function cuts the give table into its appropriate bins
  scores_x1 <- cut(as.matrix(x1),bins,labels = c('0','1','10','1E2','1E3','1E4','1E5','1E6','1E7','1E8','1E9','1E10')) 
  
  #transform function convert the tables into a column format: easy for visualization 
  Table_x1<-transform(table(scores_x1)) #contains 2 columns: "scores_x1", "Freq"
  
  #Repeating the same steps for x2
  scores_x2 <- cut(as.matrix(x2),bins,labels = c('0','1','10','1E2','1E3','1E4','1E5','1E6','1E7','1E8','1E9','1E10'))
  Table_x2 <- transform(table(scores_x2))
  
  #Getting the names of x1 and x2
  arg1 <- deparse(substitute(x1))
  arg2 <- deparse(substitute(x2))
  
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
