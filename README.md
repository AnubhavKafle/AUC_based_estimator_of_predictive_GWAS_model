# Simulate GWAS data 

### Python version of the G-WIZ tool for simulating GWAS data for different types of genetic models:
(a) Dominant
(b) Additive

Original paper:
Patron J, Serra-Cayuela A, Han B, Li C, Wishart DS (2019) Assessing the performance of genome-wide association studies for predicting disease risk. PLoS ONE 14(12): e0220215. https://doi.org/10.1371/journal.pone.0220215

URL: [paper link] {https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0220215#references}
Code: [R source code] {https://github.com/jonaspatronjp/GWIZ-Rscript/}

#to run my tool
`python SNP_simulation_model.py --PRS_file <summary GWAS file> --outpath AUC_based_estimator_of_predictive_GWAS_model --Ncases 50000 --Ncontrol 50000`

### Required
(1) Python v3.5+
(2) NumPy v1.1+
(3) Pandas v0.22+
