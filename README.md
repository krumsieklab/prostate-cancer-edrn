# EDRN Prostate Cancer Metabolomics 

This repository contains the supporting code for the EDRN Prostate Cancer Metabolomics data described in Benedetti et al., Plasma Metabolomics Profiling of 580 Patients from the Weill Cornell EDRN Prostate Cancer Cohort. 

The script _EDRN_Preprocessing.R_ allows to load and preprocess the metabolomics data described in the paper, while the script _EDRN_SamplePipeline.R_ illustrates how to use the maplet R package to run modular differential analysis. As an example, we analyzed the association of metabolite levels with patients' age.

Both script will generate an html report illustrating all the analysis steps and results.

The necessary data will be automatically downloaded from the Figshare repository [link] using the function contained in _HelperFunction.R_.
