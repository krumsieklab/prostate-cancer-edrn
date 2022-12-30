#### Initialize ----

# set working director to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load libraries
library(dplyr)
library(magrittr)
library(maplet)

# source helper function
source("HelperFunction.R")


#### Define global parameters ----

# max missingness for quotient normalization reference sample calculation
max_miss_norm <- 0.2
# max missingness for filtering
max_miss_filter <- 0.5
# significance threshold
alpha <- 0.05

#### Download Data Files ----

# download metabolomics data from figshare
# data files will be saved in the working directory
# If file already exists in working directory, will use local copy.
file_data <- "EDRN_MetabolomicsData.xlsx"
load.web.file(
  url="https://figshare.com/ndownloader/files/38682632?private_link=90221d74cdf8ebad84c5",
  md5sum = "4a5b92b24615a6a0499804b4842600b9",
  outfile = file_data
)

# download clinical data from figshare
# data files will be saved in the working directory
# If file already exists in working directory, will use local copy.
file_clin <- "EDRN_ClinicalData.xlsx"
load.web.file(
  url="https://figshare.com/ndownloader/files/38677247?private_link=90221d74cdf8ebad84c5",
  md5sum = "47c9db61e17f8ddd5ce3ef28d7f2f68c",
  outfile = file_clin
)


#### Lod Data ----

D <- 
  # load unprocessed metabolomics data
  mt_load_xls(file=file_data, sheet="Data", samples_in_row=T, id_col="SAMPLE_ID") %>% 
  mt_anno_xls(file=file_data, sheet="SampleAnnotations", anno_type="samples", anno_id_col ="SAMPLE_ID", data_id_col = "SAMPLE_ID") %>% 
  mt_anno_xls(file=file_data, sheet="MetaboliteAnnotations", anno_type="features", anno_id_col="BIOCHEMICAL", data_id_col = "name") %>%
  # load clinical data
  mt_anno_xls(file=file_clin, sheet="ClinicalAnnotations", anno_type="samples", anno_id_col = "SAMPLE_ID", data_id_col = "SAMPLE_ID") %>%
  # set variable to factor
  mt_anno_apply(anno_type = "samples", col_name = "Diagnosis", fun = as.factor)


#### Preprocess Data ----

D %<>%
  # header for html report
  mt_reporting_heading(heading = "Preprocessing", lvl = 1) %>%
  # plot metabolite missingness
  mt_plots_missingness(feat_max=max_miss_filter) %>%
  # plot sample boxplot
  mt_plots_sample_boxplot(color=Diagnosis, title = "Original", plot_logged = T) %>%
  # normalize with probabilistic quotient method
  mt_pre_norm_quot(feat_max = max_miss_norm) %>% 
  # plot dilution factor
  mt_plots_dilution_factor(in_col="Diagnosis") %>%
  # plot sample boxplot after normalization
  mt_plots_sample_boxplot(color=Diagnosis,title = "After normalization", plot_logged = T) %>%
  # log transform
  mt_pre_trans_log() %>%

  # header for html report
  mt_reporting_heading(heading = "Fisher's Missingness Analysis, Diagnosis", lvl = 1) %>%
  # Fisher's test on missing metabolites vs Diagnosis
  mt_stats_univ_missingness(in_col="Diagnosis", stat_name="missDiagnosis") %>%
  # correct test p-values for multiple tests
  mt_post_multtest(stat_name="missDiagnosis", method="BH") %>%
  # report summary statistics in html
  mt_reporting_stats(stat_name = "missDiagnosis", stat_filter = p.adj < alpha) %>%

  # header for html report
  mt_reporting_heading(heading = "Filtering", lvl = 1) %>%
  # filter metabolites with more than max_miss_filter missing values
  mt_pre_filter_missingness(feat_max=max_miss_filter) %>%
  # plot metabolite missingness
  mt_plots_missingness(feat_max=max_miss_filter) %>%
  # knn-imputation, takes several minutes
  mt_pre_impute_knn() %>%
  # plot sample boxplots
  mt_plots_sample_boxplot(color=Diagnosis, title = "After imputation", plot_logged = F)
  

#### Global Statistics ----

D %<>%
  # header for html report
  mt_reporting_heading(heading="Global Statistics", lvl=1) %>% 
  # plot PCA
  mt_plots_pca(scale_data = T, title = "scaled PCA - Diagnosis", color=Diagnosis, size=2.5, ggadd=ggplot2::scale_size_identity()) %>%
  mt_plots_pca(scale_data = T, title = "scaled PCA - Age", color=Age, size=2.5, ggadd=ggplot2::scale_size_identity()) %>%
  # plot UMAP
  mt_plots_umap(scale_data = T, title = "scaled UMAP - Diagnosis", color=Diagnosis, size=2.5, ggadd=ggplot2::scale_size_identity()) %>%
  mt_plots_umap(scale_data = T, title = "scaled UMAP - Age", color=Age, size=2.5, ggadd=ggplot2::scale_size_identity()) %>%
  # plot dataset heatmap
  mt_plots_heatmap(scale_data = T, fontsize = 5, 
                   annotation_col = c("Diagnosis","Age"), 
                   annotation_row = c("SUPER_PATHWAY"), 
                   show_colnames = F, show_rownames=F,
                   clustering_method = "ward.D2", cutree_rows = 5, cutree_cols = 5)


#### Write Preprocessed Data and Annotations to File ----
  
D %>% 
  mt_write_se_xls(file="EDRN_MetabolomicsData_Preprocessed.xlsx",
                      sheet_names = c("PreprocessedData","MetaboliteAnnotations","ClinicalAnnotations"))


#### Save html report ----

# takes several minutes
D %<>%
  mt_reporting_html(file = "EDRN_Preprocessing.html")


#### Access Data ----

# extract data matrix
dt <- D %>% assay
# extract sample annotations
cd <- D %>% colData %>% as.data.frame
# extract metabolite annotations
rd <- D %>% rowData %>% as.data.frame
