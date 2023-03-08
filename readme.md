<font size =3> To be edited</font>:

This repository contains code and example data for LC-MS/MS and amplicon sequencing data integration and analysis of the CCE1706 dataset.
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Functional-Metabolomics-Lab/CCE_Data-Analysis/blob/main/)

List of Notebooks provided in the GitRepo and its use:
1. ASV processing from raw data
2. Getting the elemental composition from SIRIUS CANOPUS info Van-Krevelen Plots
3. Data cleanup of feature table (combine SIRIUS, Canopus info) and PCoA & PermANOVA
4. After PCoA, we went Stacked Bar Plots (with Phyla level binning for Sequencing data; and SuperClass level binning for Metabolites data)
5. Sunburst visualisation (in Python Notebook) - until Genus level for sequencing data
6. Correlation Analysis script: (Focus on explanation of FDR, how you do it and justification for the method)
7. To list the important features, we tried RF and XGBoost (But we should focus on surface samples with different cycles)

## FBMN links:
[Analogs OFF](https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task=be9f2757d99148cc952bb5237096c7fd),
[Analogs ON](https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task=9d10e569e4254990b26b655b45f6eba7#)
