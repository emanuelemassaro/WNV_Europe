Spatial role of land cover on West Nile virus disease in Europe

This repository contains the analysis code, data, and figures associated with the manuscript:

Spatial role of land cover on West Nile virus disease in Europe

Authors:
Nicola Rizzetti, Angela Fanelli, Wojciech Szewczyk*, Alessandro Cescatti, Juan Carlos Ciscar, Grégoire Dubois, Dolores Ibarreta, Jordi Figuerola, and Emanuele Massaro

Affiliations:
a European Commission, Joint Research Centre (JRC), Via E. Fermi 2749, 21027 Ispra (VA), Italy
b European Commission, Joint Research Centre (JRC), C. Inca Garcilaso 3, 41092 Sevilla, Spain
c Estación Biológica de Doñana, Consejo Superior de Investigaciones Científicas (CSIC), Avda. Américo Vespucio 26, 41092 Sevilla, Spain
d Centro de Investigación Biomédica en Red Epidemiología y Salud Pública (CIBERESP), Avda. Monforte de Lemos 3–5, 28029 Madrid, Spain

Corresponding author:
Wojciech Szewczyk
📧 Wojciech.SZEWCZYK@ec.europa.eu

📌 Overview

This repository contains all materials needed to reproduce the spatial statistical analyses and figures for the study “Spatial role of land cover on West Nile virus disease in Europe”.

The analysis explores the relationship between land cover characteristics and West Nile virus (WNV) incidence across Europe, using spatial econometric models, Geographically Weighted Regression (GWR), and a suite of sensitivity analyses.

The structure of the repository follows principles of scientific reproducibility.

📁 Repository Structure
├── Codes/
│   ├── modelsALL.ipynb
│   └── sensitivity.ipynb
│
├── Data/
│   ├── input_data/          # All spatial and tabular data used in the analysis
│   ├── processed_data/      # Aggregated or transformed datasets
│   └── model_outputs/       # Intermediate results and outputs for figures/tables
│
└── Figures/
    ├── Main/                # Figures included in the main manuscript
    └── Supplementary/       # Figures used in SI, diagnostics, and sensitivity analyses

📚 Contents
1. Codes/

This directory contains the Jupyter notebooks used to run the full analysis pipeline.

modelsALL.ipynb

Main analysis for the manuscript

Preprocessing of covariates

Spatial weight matrix construction

OLS and Moran’s I diagnostics

Geographically Weighted Regression (GWR)

Model outputs and generation of figures for the main paper

sensitivity.ipynb

Sensitivity analyses of:

Distance-band threshold selection (100–200 km)

Spatial autocorrelation patterns

Multicollinearity diagnostics (correlation matrix, VIF)

Robustness of covariate selection

Figures and tables used in the Supplementary Information

2. Data/

Contains all the data needed to reproduce the results.

Possible subfolders (depending on your structure):

input_data/
Raw spatial layers, land cover information, epidemiological data, and socio-economic variables.

processed_data/
Harmonized NUTS-level datasets, standardized covariates, and intermediate files created during the analysis.

model_outputs/
Outputs from OLS, Moran’s I, GWR fits, and sensitivity analyses.

⚠️ Some datasets may not be publicly available due to licensing or confidentiality. When this is the case, placeholders and instructions are provided.

3. Figures/

All figures generated during the analysis.

Main/

Figures included in the primary manuscript such as:

Maps of WNV incidence

Spatial patterns of land cover

GWR coefficient surfaces

Model diagnostics

Supplementary/

Figures for the SI:

Distance-band sensitivity plots

Correlation heatmaps

VIF and multicollinearity diagnostics

Residual Moran’s I across distance thresholds

⚙️ Reproducibility and Requirements
Software

Python ≥ 3.9

Jupyter Notebook

Recommended environment: conda or mamba

Key Python Libraries
geopandas  
pandas  
numpy  
matplotlib  
seaborn  
libpysal  
esda  
mgwr  
statsmodels  
scikit-learn


To install the full environment:

pip install -r requirements.txt


(or provide your full list if available)

🚀 Running the Analysis

Clone the repository:

git clone https://github.com/<your-repo-name>.git
cd <your-repo-name>


Ensure all datasets are placed in the /Data/ folder as indicated.

Open and run the notebooks in the /Codes/ directory in the following order:

modelsALL.ipynb

sensitivity.ipynb

Figures will be automatically exported to the /Figures/ directory.

📄 Citation

If you use this code or data, please cite:

Riccetti, N., Fanelli, A., Szewczyk, W.*, Cescatti, A., Ciscar, J.C., Dubois, G., Ibarreta, D., Figuerola, J., & Massaro, E.
Spatial role of land cover on West Nile virus disease in Europe.
European Commission, Joint Research Centre (JRC), 2024.

📝 License

Specify your license, e.g.:

MIT License

Creative Commons Attribution 4.0

“Available for research purposes only”

🤝 Contact

For questions related to the analysis, methodology, or data:

📧 Wojciech Szewczyk
Wojciech.SZEWCZYK@ec.europa.eu

For repository or code-related issues:

📧 Emanuele Massaro
(Insert preferred academic or personal email if desired)

✔️ Done!

Let me know if you'd like:

A shorter README for GitHub

A more technical one for reproducibility packages

A badge layout (Zenodo DOI, binder, etc.)
