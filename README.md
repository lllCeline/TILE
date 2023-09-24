# TILE
## 1.Introduction
The **TILE** package provides tools to infer the fractions of a priori known cell subtypes present in a sample representing a mixture of such cell-types. Inference proceeds via one of 4 methods (linear model, ssGSEA, simple average, geometric mean).

Besides, for comparison of two patient groups, we provide a variety of statistical methods and visualization methods.

For now, the package contains **7 priori known signatures** and **2 public databases**.

## 2.Installation
The dependencies includs GSVA, dplyr, ggplot2, ggpubr, ComplexHeatmap, grid, limma, reshape2, rstatix
```R
devtools::install("TILE")
library R package

library("TILE")
```
## 3.Manual
### 3.1 Cell type fractions estimtation
#### 3.1.1 CYT signature
Cytolytic activity is associated with counter-regulatory immune responses and improved prognosis, the CYT (cytolytic activity) marker is estimated using the expression of granzyme A and perforin according to Rooney et al.,2015
```R
result <- CYT_calculation(demo_TPM)

#cell type
CYT_value
PRF1
GZMA
```
#### 3.1.2 TILs marker
Previous study found tumor infiltrating lymphocytes (TILs) can inform immune oncology research and the choice of immunotherapy for individual patients.

In addition, cell population relative to total TILs cell type enrichment scores has prognostic relevance, which can be interpreted as measuring the abundance or depletion of each cell population relative to total TILs.

Total TILs score and cell type enrichment scores are calculated according to Danaher et al.,2017.
```R
result <- Co_expression(demo_TPM)
TILs_score <- data.frame(result[[1]],check.names = FALSE)
cell_enrich <- data.frame(result[[2]],check.names = FALSE)

#cell type  
B-cells  
CD45  
CD8 T cells  
Cytotoxic cells  
DC  
Exhausted CD8  
Macrophages  
Mast cells  
NK CD56dim cells  
NK cells  
Neutrophils  
T-cells  
Th1 cells  
Treg  
TILs_Score  
```
#### 3.1.3 TILs_28_subpopulations
Pornpimol et al.,2017 selected 782 genes which are representative of specific immune cell subpopulations.

We then used gene set enrichment analysis (Barbie, David A et al.2009) with this gene set to decompose cellular profiles from RNA sequencing data from bulk tissue for individual samples.
```R
#input matrix: raw reads count
result <- ssgsea(demo_counts,normalization = FALSE,db = "TILs_28_subpopulations")
#input matrix: normalized expression values
result <- ssgsea(demo_TPM,normalization = TRUE,db = "TILs_28_subpopulations")

#cell type
Activated_B_cell   
Activated_CD4_T_cell  
Activated_CD8_T_cell  
Activated_dendritic_cell  
CD56bright_natural_killer_cell  
CD56dim_natural_killer_cell  
Central_memory_CD4_T_cell  
Central_memory_CD8_T_cell  
Effector_memeory_CD4_T_cell  
Effector_memeory_CD8_T_cell  
Eosinophil  
Gamma_delta_T_cell  
Immature_B_cell  
Immature_dendritic_cell  
MDSC  
Macrophage  
Mast_cell  
Memory_B_cell  
Monocyte  
Natural_killer_T_cell  
Natural_killer_cell  
Neutrophil  
Plasmacytoid_dendritic_cell  
Regulatory_T_cell  
T_follicular_helper_cell  
Type_17_T_helper_cell  
Type_1_T_helper_cell  
Type_2_T_helper_cell 
```
#### 1.3.4 IPRES marker
Innately resistant tumors display a transcriptional signature (referred to as the IPRES or Innate anti-PD-1 Resistance) indicating concurrent upexpression of genes involved in the regulation of mesenchymal transition, cell adhesion, ECM remodeling, angiogenesis and wound-healing (Hugo et al.,2016).

Previous findings suggest that attenuating the biological processes that underlie IPRES may improve anti-PD1 response in melanoma and other cancer types.
```R
#input matrix: raw reads count
result <- IPRES(demo_counts,normalization = FALSE)
#input matrix: normalized expression values
result <- IPRES(demo_TPM)

#cell type  
ANASTASSIOU_MULTICANCER_INVASIVENESS_SIGNATURE  
CHARAFE_BREAST_CANCER_BASAL_VS_MESENCHYMAL_DN  
EP_BLOOD_VESS_DEVEL_DN_IN_R  
HARRIS_HYPOXIA  
JAEGER_METASTASIS_UP  
JEON_SMAD6_TARGETS_DN  
KARAKAS_TGFB1_SIGNALING  
LEF1_UP.V1_UP  
LIEN_BREAST_CARCINOMA_METAPLASTIC  
LU_TUMOR_ANGIOGENESIS_UP  
LU_TUMOR_ENDOTHELIAL_MARKERS_UP  
LU_TUMOR_VASCULATURE_UP  
MAHADEVAN_GIST_MORPHOLOGICAL_SWITCH  
MAINA_VHL_TARGETS_DN  
MAPKi_INDUCED_ANGIOGENESIS  
MAPKi_INDUCED_EMT  
MISHRA_CARCINOMA_ASSOCIATED_FIBROBLAST_UP  
MS_RESP_TO_HYPOXIA_UP_IN_MAPKi_aPDL1_NR  
MS_RESP_TO_WOUNDING_UP_IN_MAPKi_aPDL1_NR  
POOLA_INVASIVE_BREAST_CANCER_UP  
POST_OP_WOUNDHEALING  
ROY_WOUND_BLOOD_VESSEL_UP  
VECCHI_GASTRIC_CANCER_ADVANCED_VS_EARLY_UP  
WESTON_VEGFA_TARGETS_12HR  
WESTON_VEGFA_TARGETS_6HR  
YE_METASTATIC_LIVER_CANCER  
```
#### 1.3.5 GEP marker
Mark et al.,2017 identified a biomarker T-cell–inflamed gene-expression profile (GEP), which was tumor type–independent dimensions of the tumor microenvironment relevant to predicting clinical outcome for agents targeting the PD-1/PD-L1 signaling pathway.
```R
result <- GEP(exp_matrix = demo_TPM)

#cell type  
GEP_score  
```
#### 1.3.6 integration signature
we integrate tumor microenvironment-related signatures, evaluated the abundance of immune cell populations in the tumor microenvironment according to Pornpimol et al.,2017 and Ott et al.,2019.
```R
result <- sig_calculate(exp_matrix = demo_TPM,db = "integration_signature")

#cell type  
Cytolytic-activity  
IFNg-signature  
Immunocostimulators  
Immunoinhibitors  
MHC-Class-I  
MHC-Class-II  
Non-Class  
adhesion-molecules  
chemokines 
```
#### 1.3.7 TLS signature
Tertiary lymphoid structures (TLSs) provide an immunological antineoplastic effect, tls_sig signature can be used to calculate the score of Tertiary lymphoid structures according to Tokunaga et al.,2020 and Rita et al.,2020.
```R
result <- sig_calculate(exp_matrix = demo_TPM,db = "tls_sig")

#cell type  
TLS_12gene  
TLS_9gene
```
#### 1.3.8 public databases
msigdb : https://www.gsea-msigdb.org/gsea/msigdb
```R
#db : msigdb_C2, msigdb_C5, msigdb_C6, msigdb_C7, msigdb_H
#input matrix: raw reads count
result <- ssgsea(demo_counts,normalization = FALSE,db = "msigdb_H")
#input matrix: normalized expression values
result <- ssgsea(demo_TPM,normalization = TRUE,db = "msigdb_H")
```
cancerSEA : http://biocc.hrbmu.edu.cn/CancerSEA/
```R
#input matrix: raw reads count
result <- ssgsea(demo_counts, normalization = FALSE, db = "cancerSEA")
#input matrix: normalized expression values
result <- ssgsea(demo_TPM,normalization = TRUE,db = "cancerSEA")
```
### 3.2 Statistical methods
#### 3.2.1 diffES
For comparison of two patient groups, different immune cell types were analyzed by the limma package.
```R
result <- diffES(score_matrix = cell_enrich,groupinfo = demo_groupinfo)
```
#### 3.2.2 t test
For comparison of two groups,t test can be used when your data are independent and normally distributed.

Parameter paired is setted as "True", when compare the means of two samples when each observation in one sample can be paired with an observation in the other sample.
```R
result <- t_test(score_matrix = cell_enrich,groupinfo = demo_groupinfo,paired = F)
#paired sample
result <- t_test(score_matrix = cell_enrich,groupinfo = simulated_groupinfo_paired,paired = T)
```
#### 3.2.3 wilcoxon test
For comparison of two groups,non-parametric two-sided Wilcoxon-rank sum test was used.

Parameter paired is setted as "True", when compare the means of two samples when each observation in one sample can be paired with an observation in the other sample.
```R
result <- wilcox(score_matrix = cell_enrich,groupinfo = demo_groupinfo,paired = F)
#paired sample
result <- wilcox(score_matrix = cell_enrich,groupinfo = simulated_groupinfo_paired,paired = T)
```
### 3.3 Visualization methods  
```R
###box plots for immnume signatures score with clinical features
p <- marker_boxplot(score_matrix = cell_enrich,groupinfo = demo_groupinfo,fontsize = 10,paired = FALSE)
print(p)
#paired sample
p <- marker_boxplot(score_matrix = cell_enrich,groupinfo = simulated_groupinfo_paired,paired = T)
print(p)
###specify cell type variables for faceting the plot into multiple panels
p <- boxplot_diff(score_matrix = cell_enrich,groupinfo = demo_groupinfo,paired = F)
print(p)
#paired sample
p <- boxplot_diff(score_matrix = cell_enrich,groupinfo = simulated_groupinfo_paired,paired = T)
print(p)
###supervised analyses of cell type score between groups
p <- marker_heatmap(score_matrix = cell_enrich,scaled = "row",groupinfo = groupinfo,title_size = 8)
print(p)
```
### 3.4 Other tools  
Normalized gene expression data are needed as input matirx for 5/8 signatures.
This function provides expression data normalization methods.
```R
#gene_len can be one of hg19_gene_symbol_length, hg38_gene_symbol_length
#normalization method can be FPKM or TPM
result <- counts_normlization(demo_counts,method = "FPKM",gene_len = "hg38_gene_symbol_length")
```
