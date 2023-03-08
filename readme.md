# _To be edited_:

This repository contains code and example data for LC-MS/MS and amplicon sequencing data integration and analysis of the CCE1706 dataset.
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Functional-Metabolomics-Lab/CCE_Data-Analysis/blob/main/)

---

**List of Notebooks provided in this repository and its use:**
1. ASV_raw_data_processing.ipynb <br>
     - The ASV 16S and 18SV9 files are used here. The sample names are renamed similar to that of the mzxml file names. Also, the taxonomic information for each feature is combined with the OTU table.
2. Elemental_composition.ipynb <br>
     - The summary files (.tsv format) from SIRIUS and CANOPUS outputs are used to get the elemental composition and some elemental ratios for each feature. Also van-Krevelen plots can be made.
3. Data cleanup of feature table (combine SIRIUS, Canopus info) and PCoA & PermANOVA<br>
     - _to be updated_
5. After PCoA, we went Stacked Bar Plots (with Phyla level binning for Sequencing data; and SuperClass level binning for Metabolites data)<br>
     - _to be updated_
6. Sunburst visualisation (in Python Notebook) - until Genus level for sequencing data <br>
     - _to be updated_
7. Correlation Analysis script: (Focus on explanation of FDR, how you do it and justification for the method)<br>
     - _to be updated_
8. To list the important features, we tried RF and XGBoost (But we should focus on surface samples with different cycles) <br>
     - _to be updated_

## FBMN links:
[Analogs OFF](https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task=be9f2757d99148cc952bb5237096c7fd),
[Analogs ON](https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task=9d10e569e4254990b26b655b45f6eba7#)
