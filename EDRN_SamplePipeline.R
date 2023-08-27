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

# significance threshold
alpha <- 0.05


#### Download Data Files ----

file_data <- "EDRN_MetabolomicsData_2023_08_27.xlsx"
file_clin <- "EDRN_ClinicalData.xlsx"

# # download metabolomics data from figshare
# # data files will be saved in the working directory
# # If file already exists in working directory, will use local copy.
# load.web.file(
#   url="https://figshare.com/ndownloader/files/38682632?private_link=90221d74cdf8ebad84c5",
#   md5sum = "4a5b92b24615a6a0499804b4842600b9",
#   outfile = file_data
# )
# 
# # download clinical data from figshare
# # data files will be saved in the working directory
# # If file already exists in working directory, will use local copy.
# load.web.file(
#   url="https://figshare.com/ndownloader/files/38677247?private_link=90221d74cdf8ebad84c5",
#   md5sum = "47c9db61e17f8ddd5ce3ef28d7f2f68c",
#   outfile = file_clin
# )


#### Lod Data ----

D <- 
  # load preprocessed metabolomics data
  mt_load_xls(file=file_data, sheet="PreprocessedData", samples_in_row=T, id_col="CLIENT_SAMPLE_ID") %>% 
  mt_anno_xls(file=file_data, sheet="SampleAnnotations", anno_type="samples", anno_id_col ="CLIENT_SAMPLE_ID", data_id_col = "CLIENT_SAMPLE_ID") %>% 
  mt_anno_xls(file=file_data, sheet="MetaboliteAnnotations", anno_type="features", anno_id_col="CHEMICAL_NAME", data_id_col = "name") %>%
  # load clinical data
  mt_anno_xls(file=file_clin, sheet="ClinicalAnnotations", anno_type="samples", anno_id_col = "SAMPLE_ID", data_id_col = "CLIENT_SAMPLE_ID") %>%
  # set variable to factor
  mt_anno_apply(anno_type = "samples", col_name = "Diagnosis", fun = as.factor)

#### Sample Pipeline ----

# Compute association of metabolite levels with Age using a linear model
D %<>%
  # header for html report
  mt_reporting_heading(heading = "Age association", lvl=1) %>%
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
                   feat_filter = p.adj < 1E-10,
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
                             color_col = "SUPER_PATHWAY")


#### Save html report ----

# takes several minutes
D %<>%
  mt_reporting_html(file = "EDRN_SamplePipeline.html")


#### Access Data ----

# extract data matrix
dt <- D %>% assay
# extract sample annotations
cd <- D %>% colData %>% as.data.frame
# extract metabolite annotations
rd <- D %>% rowData %>% as.data.frame