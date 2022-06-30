############## Indexing ASV files ###############


#### First, 16S data ####


## Call libraries needed for R
library(dplyr)
library(tidyr)

##  1) Set working directory and pull files into R

setwd("~/Graduate School/Projects/CCE-LTER Cruise/P1706/P1706_ASV_index")

table16S <- read.csv('16S_ASV_table_noplastid_newID (1).csv')
ASVtax16S <- read.csv('16S_ASV_taxonomy.csv')

##  2) Next, index files based on Feature.ID column

# Rename column X to Feature.ID to match other excel file
table16S <- rename(table16S, 'Feature.ID' = 'X')

# Remove the A_ header in all the name in the 16S table
table16S$Feature.ID <- gsub('A_', '', as.character(table16S$Feature.ID))

# Merge two datafiles by Feature.ID
ASV16S <- merge(table16S, ASVtax16S, by = 'Feature.ID')

#   3) Next, we must index the Taxon column containing multiple assignments into individual columns 

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

## Write dataframe into .csv
#Use this file when merging ASV16S data with mzmine features.
write.csv(ASV16S, file = 'P1706_ASV16S.csv')


#### 4) Rename ASV file column headers to match with mzmine column headers (sample names)

# I first had to create a .csv file that contained names of PPL samples side by side to names of genetic samples
# .csv file is 'joint names.csv'

# Rename ASV tables with mass spec sample ID. Big problem. In some samples, 1 ASV, had to manually curate
jointnames <- na.omit(read.csv('joint names.csv'))
jointnames <- jointnames %>% na_if("") %>% na.omit



# Repeat columns that have duplicate PPL samples (this way both duplicates will have ASV information tied to them)

# Merge ASV16S with mzmine df
match(jointnames[, "name1"], names(ASV16S))
names(ASV16S)[match(jointnames[, "name1"], names(ASV16S))] = jointnames[, "name4"]

colnames(ASV16S) <- sub('_MSMS', '.mzxml', as.character(colnames(ASV16S)))

write.csv(ASV16S, file = 'ASV16S.csv')

########### Aggregating into Taxonomy

# #Phylum
# 
# ASV16S_Phylum <- select(ASV16S, -c(1, 180, 182:187))
# 
# ASV16S_Phylum <- ASV16S_Phylum %>%
#   group_by(Phylum16S) %>%
#   summarize(across(starts_with('CCE'), sum))
# 
# #Class
# 
# ASV16S_Class <- select(ASV16S, -c(1, 180:181, 183:187))
# 
# ASV16S_Class <- ASV16S_Class %>%
#   group_by(Class16S) %>%
#   summarize(across(starts_with('CCE'), sum))

#Order

ASV16S_Order <- select(ASV16S, -c(1, 180:182, 184:187))

ASV16S_Order <- ASV16S_Order %>%
  group_by(Order16S) %>%
  summarize(across(starts_with('CCE'), sum))

write.csv(ASV16S_Order, file = 'ASV16S_Order.csv')

#Family

ASV16S_Family <- select(ASV16S, -c(1, 180:183, 185:187))

ASV16S_Family <- ASV16S_Family %>%
  group_by(Family16S) %>%
  summarize(across(starts_with('CCE'), sum))

write.csv(ASV16S_Family, file = 'ASV16S_Family.csv')

#Genus

ASV16S_Genus <- select(ASV16S, -c(1, 180:184, 186:187))

ASV16S_Genus <- ASV16S_Genus %>%
  group_by(Genus16S) %>%
  summarize(across(starts_with('CCE'), sum))

write.csv(ASV16S_Genus, file = 'ASV16S_Genus.csv')

# #Species
# 
# ASV16S_Species <- select(ASV16S, -c(1, 180:185, 187))
# 
# ASV16S_Species <- ASV16S_Species %>%
#   group_by(Species16S) %>%
#   summarize(across(starts_with('CCE'), sum))






################################# 18SV9 data ###########################################






###   1) Set working directory and pull files into R

setwd("~/Graduate School/Projects/CCE-LTER Cruise/P1706/P1706_ASV_index")
table18SV9 <- read.csv('18SV9_ASV_table.csv')
ASVtax18SV9 <- read.csv('18SV9_ASV_tax_table (1).csv')

###   2) Next, index files based on Feature.ID column

#Rename column X to Feature.ID to match other excel file
table18SV9 <- rename(table18SV9, 'Feature.ID' = 'X')

#merge two datafiles by Feature.ID
ASV18SV9 <- merge(table18SV9, ASVtax18SV9, by = 'Feature.ID')

###   3) Next, we must index the Taxon column containing multiple assignments into individual columns 

##  We can change the name of these headers (I didn't want to lose any information yet)
ASV18SV9 <- separate(data = ASV18SV9, col = Taxon, into = c('Domain18SV9', 'Kingdom18SV9', 
                                                            'Phylum18SV9', 'Class18SV9', 
                                                            'Order18SV9', 'Family18SV9',
                                                            'Genus18SV9', 'Species18SV9'), sep = ";")


## Write dataframe into .csv
#Use this file when merging ASV16S data with mzmine features.
 write.csv(ASV18SV9, file = 'P1706_ASV18SV9.csv')
 
 
 # Rename ASV tables with mass spec sample ID. Big problem. In some samples, 1 ASV, had to manually curate
 jointnames <- na.omit(read.csv('joint names.csv'))
 jointnames <- jointnames %>% na_if("") %>% na.omit
 
# Creating Jointnames3 to use to index 18S data
 jointnames3 <- jointnames[!duplicated(jointnames$name2),]
 names(ASV18SV9)[match(jointnames3[, "name2"], names(ASV18SV9))] = jointnames3[, "name1"]


 #This following step commented out because it is not necessary.
 #Giving repeated values in jointnames unique values.
 #jointnames2$name2 <- make.names(jointnames2$name2, unique=T)
 
 
 # Give ASV18SV9 similar file name as metadata. 
 match(jointnames[, "name1"], names(ASV18SV9))
 names(ASV18SV9)[match(jointnames[, "name1"], names(ASV18SV9))] = jointnames[, "name4"]
 
 colnames(ASV18SV9) <- sub('_MSMS', '.mzxml', as.character(colnames(ASV18SV9)))

 write.csv(ASV18SV9, file = 'ASV18SV9.csv')

 ########### Aggregating into Taxonomy

# #Kingdom
# 
# ASV18SV9_Kingdom <- select(ASV18SV9, -c(1, 180:181, 183:189))
# 
# ASV18SV9_Kingdom <- ASV18SV9_Kingdom %>%
#   group_by(Kingdom18SV9) %>%
#   summarize(across(starts_with('CCE'), sum))
# 
# #Phylum
# 
# ASV18SV9_Phylum <- select(ASV18SV9, -c(1, 180:182, 184:189))
# 
# ASV18SV9_Phylum <- ASV18SV9_Phylum %>%
#   group_by(Phylum18SV9) %>%
#   summarize(across(starts_with('CCE'), sum))
# 
# #Class
# 
# ASV18SV9_Class <- select(ASV18SV9, -c(1, 180:183, 185:189))
# 
# ASV18SV9_Class <- ASV18SV9_Class %>%
#   group_by(Class18SV9) %>%
#   summarize(across(starts_with('CCE'), sum))

#Order

ASV18SV9_Order <- select(ASV18SV9, -c(1, 180:184, 186:189))

ASV18SV9_Order <- ASV18SV9_Order %>%
  group_by(Order18SV9) %>%
  summarize(across(starts_with('CCE'), sum))

write.csv(ASV18SV9_Order, file = 'ASV18SV9_Order.csv')

#Family

ASV18SV9_Family <- select(ASV18SV9, -c(1, 180:185, 187:189))

ASV18SV9_Family <- ASV18SV9_Family %>%
  group_by(Family18SV9) %>%
  summarize(across(starts_with('CCE'), sum))

write.csv(ASV18SV9_Family, file = 'ASV18SV9_Family.csv')

#Genus

ASV18SV9_Genus <- select(ASV18SV9, -c(1, 180:186, 188:189))

ASV18SV9_Genus <- ASV18SV9_Genus %>%
  group_by(Genus18SV9) %>%
  summarize(across(starts_with('CCE'), sum))

write.csv(ASV18SV9_Genus, file = 'ASV18SV9_Genus.csv')

# #Species
# 
# ASV18SV9_Species <- select(ASV18SV9, -c(1, 180:187, 189))
# 
# ASV18SV9_Species <- ASV18SV9_Species %>%
#   group_by(Species18SV9) %>%
#   summarize(across(starts_with('CCE'), sum))