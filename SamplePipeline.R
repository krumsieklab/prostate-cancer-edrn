#### Initialize ----

# set working director to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load libraries
library(dplyr)
library(magrittr)
library(maplet)

#### Define global parameters ----

# max missingness for quotient normalization reference sample calculation
max_miss_norm <- 0.2
# max missingness for filtering
max_miss_filter <- 0.5
# significance threshold
alpha <- 0.05

#### Define Helper Functions ----

# Downloads files from the web and verifies their checksum. 
# If file already exists in working directoyy, will use local copy.
load.web.file <- function(
    url, md5sum, outfile, zipfile = F) {
  # check if local file exists
  if (file.exists(outfile)) {
    # verify checksum
    realsum <- tools::md5sum(outfile)[[1]]
    if (realsum != md5sum) stop(sprintf("Local file %s has wrong checksum: %s", outfile, realsum))
    # do not delete wrong file, it was already here before
    
  } else {
    if(zipfile){
      # download file
      temp <- tempfile()
      download.file(url,temp)
      unzip(zipfile = temp, files = outfile, exdir = ".")
    } else {
      # download file
      download.file(url, outfile)
    }
    # verify checksum
    realsum <- tools::md5sum(outfile)[[1]]
    if (realsum != md5sum) { 
      # delete wrong file
      unlink(outfile)
      stop(sprintf("Remote file %s has wrong checksum: %s", url, realsum))
    }
  }
}

#### Download Data Files ----

# download metabolomics data from figshare
# data files will be saved in the working directory
file_data <- "EDRN_MetabolomicsData.xlsx"
load.web.file(
  url="https://figshare.com/ndownloader/files/38639168?private_link=90221d74cdf8ebad84c5",
  md5sum = "97629808d6612cfe29cd84bd024ba96d",
  outfile = file_data
)

# download clinical data from figshare
# data files will be saved in the working directory
file_clin <- "EDRN_ClinicalData.xlsx"
load.web.file(
  url="https://figshare.com/ndownloader/files/38628464?private_link=90221d74cdf8ebad84c5",
  md5sum = "70324a7d65eb6a213eefdb4dca92d9ca",
  outfile = file_clin
)

#### Lod Data ----

D <- 
  # load data
  mt_load_xls(file=file_data, sheet="Data", samples_in_row=T, id_col="SAMPLE_ID") %>% 
  mt_anno_xls(file=file_data, sheet="SampleAnnotations", anno_type="samples", anno_id_col ="SAMPLE_ID", data_id_col = "SAMPLE_ID") %>% 
  mt_anno_xls(file=file_data, sheet="MetaboliteAnnotations", anno_type="features", anno_id_col="BIOCHEMICAL", data_id_col = "name") %>%
  # load clinical data
  mt_anno_xls(file=file_clin, sheet="ClinicalInfo", anno_type="samples", anno_id_col = "A1", data_id_col = "SAMPLE_ID") %>%
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
  mt_reporting_heading(heading = "Fisher's Missingness Analysis", lvl = 1) %>%
  # fisher's test on missing metabolites vs Diagnosis
  mt_stats_univ_missingness(in_col="Diagnosis", stat_name="missDiagnosis") %>%
  # correct test p-values for multiple tests
  mt_post_multtest(stat_name="missDiagnosis", method="BH") %>%
  # report summary statistics in html
  mt_reporting_stats(stat_name = "missDiagnosis", stat_filter = p.adj < alpha) %>%

  # header for html report
  mt_reporting_heading(heading = "Filtering", lvl = 1) %>%
  # plot metabolite missingness
  mt_plots_missingness(feat_max=max_miss_filter) %>%
  # filter metabolites with more than max_miss_filter missing values
  mt_pre_filter_missingness(feat_max=max_miss_filter) %>%
  {.}

#### Imputation ----

D %<>%
  # header for html report
  mt_reporting_heading(heading="Global Statistics", lvl=1) %>% 
  # knn-imputation
  mt_pre_impute_knn() 

#### Global Statistics ----

D %<>%
  # plot sample boxplots
  mt_plots_sample_boxplot(color=Diagnosis, title = "After imputation", plot_logged = T) %>%
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

#### Covariate Analysis ----

D %<>%
  # header for html report
  mt_reporting_heading(heading = "Age correction", lvl=1) %>%
  # run linear model of metabolite~Age
  mt_stats_univ_lm(formula = ~Age,
                   samp_filter = (!is.na(Age)),
                   stat_name = "Age") %>%
  # correct for multiple tests
  mt_post_multtest(stat_name = "Age", method = "BH") %>%
  # report summary statistic in report
  mt_reporting_stats(stat_name = "Age", stat_filter = p.adj < !!alpha) %>%
  # plot volcano of significant results
  mt_plots_volcano(stat_name = "Age",
                   x = statistic,
                   feat_filter = p.adj < !!alpha,
                   colour       = p.adj < !!alpha) %>%
  # plot scatter plots of significant results
  mt_plots_box_scatter(plot_type = "scatter",
                       stat_name = "Age",
                       jitter="jitter",
                       x = Age,
                       feat_filter = (p.adj < 1E-10),
                       feat_sort = p.value,
                       annotation = "{sprintf('P-value: %.2e', p.value)}\nPadj: {sprintf('%.2e', p.adj)}") %>%
  # plot pathway bar plot
  mt_plots_stats_pathway_bar(stat_list = "Age", 
                             feat_filter = p.value<!!alpha,
                             y_scale = "count", assoc_sign_col = "statistic",
                             group_col = "SUB_PATHWAY",
                             color_col = "SUPER_PATHWAY") %>%
  # correct metabolites for Age using a linear model
  mt_pre_confounding_correction(formula = ~Age)

#### Write Preprocessed Data and Annotations to File ----
  
D %>% mt_write_se_xls(file="EDRN_MetabolomicsData_Preprocessed.xlsx",
                      sheet_names = c("PreprocessedData","MetaboliteAnnotations","ClinicalAnnotations"))

#### Save html report ----

D %<>%
  mt_reporting_html(file = "EDRN_Preprocessing.html")

#### Access data ----

# extract data matrix
dt <- D %>% assay
# extract sample annotations
cd <- D %>% colData %>% as.data.frame
# extract metabolite annotations
rd <- D %>% rowData %>% as.data.frame
