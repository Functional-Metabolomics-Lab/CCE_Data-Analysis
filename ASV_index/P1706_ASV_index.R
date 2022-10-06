#' The following code organizes and combines 16S and 18S sequencing files
#' to more easily do further data analysis


########################      16S and 18SV9     #########################


## Call libraries needed for R
if (!require(dplyr)) install.packages('dplyr')
library(dplyr)
if (!require(tidyr)) install.packages('tidyr')
library(tidyr)

##  1) Set working directory and pull files into R . (Change working directory if needed)

setwd("~/Graduate School/Projects/CCE-LTER Cruise/P1706/Most up to date code/ASV_index")

## 2) Load raw files into R

# 16S raw files
table16S <- read.csv('16S_ASV_table.csv')
ASVtax16S <- read.csv('16S_ASV_taxonomy.csv')
# 18S raw files
table18SV9 <- read.csv('18SV9_ASV_table.csv')
ASVtax18SV9 <- read.csv('18SV9_ASV_taxonomy.csv')
# Jointname files for renaming ASV samples (had to manually curate to match ASV names with feature table and metadata)
jointnames <- na.omit(read.csv('joint names.csv'))
jointnames <- jointnames %>% na_if("") %>% na.omit
jointnames2 <- na.omit(read.csv('joint names2.csv'))
jointnames2 <- jointnames2 %>% na_if("") %>% na.omit


##  3) 16S Processing first

##  a) Index files based on Feature.ID column

##  Rename column X to Feature.ID to match other excel file
table16S <- rename(table16S, 'Feature.ID' = 'X')

# Remove the A_ header in all of the names in the 16S table
table16S$Feature.ID <- gsub('A_', '', as.character(table16S$Feature.ID))

# Merge two datafiles by Feature.ID
ASV16S <- merge(table16S, ASVtax16S, by = 'Feature.ID')

#   b) Next, we must index the Taxon column containing multiple assignments into individual columns 

# 0 is domain, 1 is phylum, 2 is class, 3 is order, 4 is family, 5 is genus, 6 is species
##  We can change the name of these headers (I didn't want to lose any information yet)
ASV16S <- separate(data = ASV16S, col = Taxon, into = c('Domain16S', 'Phylum16S', 
                                                     'Class16S', 'Order16S', 
                                                     'Family16S', 'Genus16S',
                                                     'Species16S'), sep = ";")

## Now remove the D_0__ header in each newly formed column (There is probably a way to clean up this code...)
ASV16S$Domain16S <- gsub('D_0__', '', as.character(ASV16S$Domain16S))
ASV16S$Phylum16S <- gsub('D_1__', '', as.character(ASV16S$Phylum16S))
ASV16S$Class16S <- gsub('D_2__', '', as.character(ASV16S$Class16S))
ASV16S$Order16S <- gsub('D_3__', '', as.character(ASV16S$Order16S))
ASV16S$Family16S <- gsub('D_4__', '', as.character(ASV16S$Family16S))
ASV16S$Genus16S <- gsub('D_5__', '', as.character(ASV16S$Genus16S))
ASV16S$Species16S <- gsub('D_6__', '', as.character(ASV16S$Species16S))

#### c) Rename ASV file column headers to match with feature table (using jointnames)

## For cycles, two samples were taken for MS analysis and only one for ASV. 
# The following code can be used if we do not want to make duplicates to match ASV files to PPL samples
# 
# match(jointnames[, "name1"], names(ASV16S))
# names(ASV16S)[match(jointnames[, "name1"], names(ASV16S))] = jointnames[, "name4"]
# 
# colnames(ASV16S) <- sub('_MSMS', '.mzxml', as.character(colnames(ASV16S)))
# 
# write.csv(ASV16S, file = 'ASV16S_noduplicates.csv')

# The following code is for making duplicates to match ASV files to PPL samples
# Repeat columns that have duplicate PPL samples (this way both duplicates will have ASV information tied to them)
ASV16S <- ASV16S %>% mutate(
  across(
    .cols = contains(c('C1', 'C2', 'C3', 'C4', 'SBB', 'B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B9'),
                     ignore.case = FALSE),
    .names = '{.col}2'
  )
)

# Merge ASV16S with mzmine df
match(jointnames2[, "name1"], names(ASV16S))
names(ASV16S)[match(jointnames2[, "name1"], names(ASV16S))] = jointnames2[, "name4"]
colnames(ASV16S) <- sub('_MSMS', '.mzxml', as.character(colnames(ASV16S)))

write.csv(ASV16S, file = 'ASV16S.csv')



##  4) 18SV9 data processing next


##  a) Next, index files based on Feature.ID column

#Rename column X to Feature.ID to match other excel file
table18SV9 <- rename(table18SV9, 'Feature.ID' = 'X')

#Merge two datafiles by Feature.ID
ASV18SV9 <- merge(table18SV9, ASVtax18SV9, by = 'Feature.ID')

###   b) Next, we must index the Taxon column containing multiple assignments into individual columns 

##  We can change the name of these headers (I didn't want to lose any information yet)
ASV18SV9 <- separate(data = ASV18SV9, col = Taxon, into = c('Domain18SV9', 'Kingdom18SV9', 
                                                            'Phylum18SV9', 'Class18SV9', 
                                                            'Order18SV9', 'Family18SV9',
                                                            'Genus18SV9', 'Species18SV9'), sep = ";")
 
#### c) Rename ASV file column headers to match with feature table (using jointnames)

## For cycles, two samples were taken for MS analysis and only one for ASV. 
# The following code can be used if we do not want to make duplicates to match ASV files to PPL samples
# Creating Jointnames3 to use to index 18S data
#  jointnames3 <- jointnames[!duplicated(jointnames$name2),]
#  names(ASV18SV9)[match(jointnames3[, "name2"], names(ASV18SV9))] = jointnames3[, "name1"]
#  
#  # Give ASV18SV9 similar file name as metadata. 
#  match(jointnames[, "name1"], names(ASV18SV9))
#  names(ASV18SV9)[match(jointnames[, "name1"], names(ASV18SV9))] = jointnames[, "name4"]
#  
#  colnames(ASV18SV9) <- sub('_MSMS', '.mzxml', as.character(colnames(ASV18SV9)))
# 
#  write.csv(ASV18SV9, file = 'ASV18SV9_noduplicates.csv')
 
#### Creating Jointnames3 to use to index 18S data

#Remove duplicated names in jointnames2 to allow merge
 jointnames3 <- jointnames2[!duplicated(jointnames2$name2),]
# Matching 18S ASV names to 16S ASV names
 names(ASV18SV9)[match(jointnames3[, "name2"], names(ASV18SV9))] = jointnames3[, "name1"]

 # The following code is for making duplicates to match ASV files to PPL samples
 # Repeat columns that have duplicate PPL samples (this way both duplicates will have ASV information tied to them)
 ASV18SV9 <- ASV18SV9 %>% mutate(
   across(
     .cols = contains(c('C1', 'C2', 'C3', 'C4', 'SBB', 'B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B9'),
                      ignore.case = FALSE),
     .names = '{.col}2'
   )
 )
 
 # Rename ASV18SV9 names with feature table / metadata names
 match(jointnames2[, "name1"], names(ASV18SV9))
 names(ASV18SV9)[match(jointnames2[, "name1"], names(ASV18SV9))] = jointnames2[, "name4"]
 colnames(ASV18SV9) <- sub('_MSMS', '.mzxml', as.character(colnames(ASV18SV9)))
 
 write.csv(ASV18SV9, file = 'ASV18SV9.csv')

