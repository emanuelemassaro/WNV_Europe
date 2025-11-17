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
Wojciech Szewczyk — Wojciech.SZEWCZYK@ec.europa.eu

Overview

This repository contains all materials needed to reproduce the spatial statistical analyses and figures for the study Spatial role of land cover on West Nile virus disease in Europe.

The analysis explores the relationship between land cover characteristics and West Nile virus (WNV) incidence across Europe, using spatial econometric models, Geographically Weighted Regression (GWR), and a suite of sensitivity analyses.

Repository Structure
├── Codes/
│   ├── modelsALL.ipynb
│   └── sensitivity.ipynb
│
├── Data/
│   ├── input_data/
│   ├── processed_data/
│   └── model_outputs/
│
└── Figures/
    ├── Main/
    └── Supplementary/

Contents
1. Codes/
modelsALL.ipynb

Main workflow for:

data preparation

construction of spatial weights

OLS and Moran’s I diagnostics

GWR modelling

generation of figures for the manuscript

sensitivity.ipynb

Sensitivity analyses:

distance-band threshold selection

multicollinearity diagnostics (correlation matrix, VIF)

robustness checks for covariate selection

supplementary figures and tables

2. Data/

Contains all datasets used in the analysis:

input_data/: raw spatial and tabular data

processed_data/: harmonized NUTS-level inputs

model_outputs/: OLS, GWR, and sensitivity results

Note: some datasets may not be openly shared due to licensing restrictions. In such cases, placeholders or instructions are provided.

3. Figures/

All figures produced by the analysis.

Main/: figures included in the manuscript

Supplementary/: figures for the SI (diagnostics, sensitivity, heatmaps, etc.)

Requirements

Python ≥ 3.9
Recommended environment: Conda / Mamba

Key libraries
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


Install with:

pip install -r requirements.txt

Running the Analysis

Clone the repository:

git clone https://github.com/<your-repo>.git
cd <your-repo>


Ensure all data are placed in Data/.

Run notebooks from the Codes/ folder:

modelsALL.ipynb

sensitivity.ipynb

Figures will be exported automatically to the Figures/ folder.

Citation

Rizzetti, N., Fanelli, A., Szewczyk, W.*, Cescatti, A., Ciscar, J.C., Dubois, G., Ibarreta, D., Figuerola, J., & Massaro, E.
Spatial role of land cover on West Nile virus disease in Europe.
European Commission, Joint Research Centre (JRC), 2024.

License

Specify your license here (MIT, CC-BY-4.0, etc.)

Contact

For scientific questions:
Wojciech Szewczyk — Wojciech.SZEWCZYK@ec.europa.eu

For code or repository issues:
Emanuele Massaro
