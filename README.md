# SARAL:Genetic Risk Prediction Using a Spatial Autoregressive Model with Adaptive Lasso 
This is a R package that implements the SARAL algorithm. SARAL is designed to build genetic risk prediction models for sequencing data. SARAL is a set-based approach. It reduces the data dimension and accumulates genetic effects within a single nucleotide variant (SNV) set. SARAL allows different SNV sets having various magnitudes and directions of effect sizes. With an adaptive lasso implemented, SARAL can shrink the effects of noise SNV sets to be zero and thus further improve prediction accuracy. 

You can use devtools::install_github("YaluWen/SARAL/SARAL.package") to install the package from R. 
