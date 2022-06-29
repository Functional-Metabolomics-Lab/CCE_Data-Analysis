# Merging all Sirius and MZmine table output

####  Notes: For normalization, df will contain molecular feature output,
####  ASV16S and ASV18SV9 should be normalized separately
####  

# Call libraries needed for R
library(dplyr)
library(tidyr)
library(readr)
library(stringr)

####  1) Set working directory and pull files into R

setwd("~/Graduate School/Projects/CCE-LTER Cruise/P1706/MZmine and Sirius/P1706_index_ralph")
ASV16S <- read.csv('ASV16S.csv')
ASV18SV9 <- read.csv('ASV18SV9.csv')
mzmine <- read.csv('CCE1706_MZmine3_GNPS_1_quant.csv')
sirius <- read.csv('formula_identifications.tsv', sep = "\t")
canopus <- read.csv('canopus_formula_summary.tsv', sep = "\t")

####  2) Cleaning up and merging files together

# Removing excess text in column headers. We can do more/less on this step depending on what we want.
colnames(mzmine) <- sub('.mzXML.Peak.area', '', as.character(colnames(mzmine)))

# Renaming certain column headers to retain different molecular formulas
sirius <- rename(sirius, 'Sirius_molecularFormula' = 'molecularFormula')
canopus <- rename(canopus, 'Canopus_molecularFormula' = 'molecularFormula')

# Merging two Canopus and Sirius dataframes by id
Siriusall <- merge(canopus, sirius, by = 'id', all = TRUE)

# Indexing featureID out of id
Siriusall <- separate(data = Siriusall, col = id, into = c(NA, NA, NA, NA, NA, 'row.ID'), sep = '_')

# Merging Siriusall and Mzmine table
df <- merge(mzmine, Siriusall, by = 'row.ID', all = TRUE)

####  3) Indexing molecular formulas into C, N, H, P, O, S, etc.

# First I will do Sirius molecular formulas
molecules <- regmatches(df$Sirius_molecularFormula, gregexpr("\\b[A-Z][a-z]*\\d*", df$Sirius_molecularFormula))
molecules <- lapply(molecules, function(a) paste0(a, ifelse(grepl("[^0-9]$", a), "1", "")))

atomcounts <- lapply(molecules, function(mol) setNames(as.integer(gsub("\\D", "", mol)), gsub("\\d", "", mol)))

atoms <- unique(unlist(sapply(atomcounts, names)))
atoms <- sapply(atoms, function(atom) sapply(atomcounts, function(a) if (atom %in% names(a)) a[atom] else 0))
rownames(atoms) <- df$Sirius_molecularFormula
atoms

Sirius_formula <- as.data.frame(atoms)
rownames(Sirius_formula) <- NULL
colnames(Sirius_formula) <- paste(colnames(Sirius_formula), "Sirius", sep = '_')

df <- cbind(df, Sirius_formula)

# Next, I will do Canopus formulas
molecules <- regmatches(df$Canopus_molecularFormula, gregexpr("\\b[A-Z][a-z]*\\d*", df$Canopus_molecularFormula))
molecules <- lapply(molecules, function(a) paste0(a, ifelse(grepl("[^0-9]$", a), "1", "")))

atomcounts <- lapply(molecules, function(mol) setNames(as.integer(gsub("\\D", "", mol)), gsub("\\d", "", mol)))

atoms <- unique(unlist(sapply(atomcounts, names)))
atoms <- sapply(atoms, function(atom) sapply(atomcounts, function(a) if (atom %in% names(a)) a[atom] else 0))
rownames(atoms) <- df$Canopus_molecularFormula
atoms

Canopus_formula <- as.data.frame(atoms)
rownames(Canopus_formula) <- NULL
colnames(Canopus_formula) <- paste(colnames(Canopus_formula), "Canopus", sep = '_')

df <- cbind(df, Canopus_formula)

# df contains only the information on molecular features


##########################################  4) Joining ASV data with Metabolite Data


# I first had to create a .csv file that contained names of PPL samples side by side to names of genetic samples
# .csv file is 'joint names.csv'

# Rename row.ID to feature.ID
df <- rename(df, 'Feature.ID' = 'row.ID')

# Rename ASV tables with mass spec sample ID. Big problem. In some samples, 1 ASV, had to manually curate
jointnames2 <- na.omit(read.csv('joint names2.csv'))
jointnames2 <- jointnames2 %>% na_if("") %>% na.omit

####### For 16S first #####

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

df2 <- merge(df, ASV16S, all = TRUE)


########  For 18S  #########


# Creating Jointnames3 to use to index 18S data
jointnames3 <- jointnames2[!duplicated(jointnames2$name2),]
names(ASV18SV9)[match(jointnames3[, "name2"], names(ASV18SV9))] = jointnames3[, "name1"]


#This following step commented out because it is not necessary.
#Giving repeated values in jointnames unique values.
#jointnames2$name2 <- make.names(jointnames2$name2, unique=T)


# Repeat columns that have duplicate PPL samples
ASV18SV9 <- ASV18SV9 %>% mutate(
  across(
    .cols = contains(c('C1', 'C2', 'C3', 'C4', 'SBB', 'B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B9'),
                     ignore.case = FALSE),
    .names = '{.col}2'
  )
)

# Merge ASV18SV9 with df2 by matching names together
match(jointnames2[, "name1"], names(ASV18SV9))
names(ASV18SV9)[match(jointnames2[, "name1"], names(ASV18SV9))] = jointnames2[, "name4"]

df3 <- merge(df2, ASV18SV9, all = TRUE)

# All data should be in dataframe now, next we can merge metadata

#####    5) Merging metadata into dataframe by transposing

# Pull metadata in and clean up column headers and file names
metadata <- read.csv('metadata_CCE.csv')
colnames(metadata) <- sub('ATTRIBUTE_', '', as.character(colnames(metadata)))
metadata$filename <- gsub('Blank_', '', as.character(metadata$filename))
metadata$filename <- gsub('.mzxml', '_MSMS', as.character(metadata$filename))

# Transpose dataframe by remembering filename and transposing all but the first column (name)
n <- metadata$filename
metadata2 <- as.data.frame(t(metadata[,-1]))
colnames(metadata2) <- n
metadata2$metadata <- factor(row.names(metadata2))

# Next, merge datafiles, metadata and feature ID condensed into column from all rows. 

wholedf <- merge(metadata2, df, all = TRUE)







#  Write dataframe into .csv
write.csv(wholedf, file = 'molecular_information.csv')
