# Simulate GWAS data for case-control settings

### Python version of the G-WIZ tool for simulating GWAS data for different types of genetic models:
1. Dominant
2. Additive

Original paper:
Patron J, Serra-Cayuela A, Han B, Li C, Wishart DS (2019) Assessing the performance of genome-wide association studies for predicting disease risk. PLoS ONE 14(12): e0220215. https://doi.org/10.1371/journal.pone.0220215

URL: [paper link] (https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0220215#references)

Code: [R source code] (https://github.com/jonaspatronjp/GWIZ-Rscript/)

# Command used to run the program

`python SNP_simulation_model.py --PRS_file <summary GWAS file> --outpath <output folder> --Ncases <N1> --Ncontrol <N2>`

Type `python SNP_simulation_model.py --help` to display all the required parameters.

### Required dependencies
- Python v3.5+
- NumPy v1.1+
- Pandas v0.22+
