# Spatial role of land cover on West Nile virus disease in Europe

This repository contains the analysis code, data, and figures associated with the manuscript:

**Spatial role of land cover on West Nile virus disease in Europe**

**Authors:**  
Nicola Riccetti, Angela Fanelli, Wojciech Szewczyk\*, Alessandro Cescatti, Juan Carlos Ciscar, Grégoire Dubois, Dolores Ibarreta, Jordi Figuerola, and Emanuele Massaro

**Affiliations:**  
a European Commission, Joint Research Centre (JRC), Via E. Fermi 2749, 21027 Ispra (VA), Italy  
b European Commission, Joint Research Centre (JRC), C. Inca Garcilaso 3, 41092 Sevilla, Spain  
c Estación Biológica de Doñana, Consejo Superior de Investigaciones Científicas (CSIC), Avda. Américo Vespucio 26, 41092 Sevilla, Spain  
d Centro de Investigación Biomédica en Red Epidemiología y Salud Pública (CIBERESP), Avda. Monforte de Lemos 3–5, 28029 Madrid, Spain

**Corresponding author:**  
*Wojciech Szewczyk* — Wojciech.SZEWCZYK@ec.europa.eu

---

## Overview

This repository contains all materials needed to reproduce the spatial statistical analyses and figures for the study *Spatial role of land cover on West Nile virus disease in Europe*.  
It includes the full codebase, processed datasets, model outputs, and figures for both the main manuscript and the Supplementary Information (SI).

---

## Repository Structure

```
├── Codes/
│   ├── modelsALL.ipynb        # Main analysis and manuscript figures
│   └── sensitivity.ipynb      # Sensitivity analyses and diagnostics
│
├── Data/
│   ├── input_data/            # Raw spatial and epidemiological data
│   ├── processed_data/        # Harmonized and aggregated datasets
│   └── model_outputs/         # Outputs from models and diagnostics
│
└── Figures/
    ├── Main/                  # Figures included in the manuscript
    └── Supplementary/         # Figures for the SI (diagnostics, heatmaps, etc.)
```
## Note on `modifyIncidenceData.ipynb` (Data Privacy)

The notebook `modifyIncidenceData.ipynb` anonymizes the original **Incidence** values from ECDC (Tessy) to comply with the data-sharing agreement, which does not allow publication of case-based data.

The notebook:
- replaces the original incidence values with **synthetic, privacy-safe values**,  
- applies **local reshuffling** between nearby rows to remove any link to the original data.

The resulting files preserve only high-level statistical patterns and are safe to publish.  
For access to the real data or more information, please contact **ECDC**.


---

## Contents

### 1. Codes/

#### `modelsALL.ipynb`
Main workflow including:
- Data preparation and standardization  
- Construction of spatial weights matrices  
- OLS regression and Moran’s I diagnostics  
- Geographically Weighted Regression (GWR)  
- Generation of main manuscript figures  

#### `sensitivity.ipynb`
Includes:
- Distance-band threshold sensitivity (100–200 km)  
- Multicollinearity diagnostics (correlation matrix, VIF)  
- Robustness tests for covariate selection  
- Supplementary figures and tables

#### `APES_clima.R`  `APES_Base 1.R`
Includes:
- R codes to aggregated the climatic, land use and infection data at NUTS 3 level

---

### 2. Data/

Includes all datasets required to run the models:

- `input_data/`: raw spatial layers and source data  
- `processed_data/`: cleaned and aligned NUTS-level datasets  
- `model_outputs/`: OLS, Moran’s I, GWR results, sensitivity results  

> Some datasets may not be publicly distributed due to licensing restrictions.  
> When this occurs, placeholder files and instructions are provided.

---

### 3. Figures/

- **Main/**: All figures used in the primary manuscript  
- **Supplementary/**: Heatmaps, residual diagnostics, sensitivity charts, and additional spatial visualizations

---

## Requirements

Python ≥ 3.9  
Recommended environment: Conda / Mamba

### Required libraries
```
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
```

Install via:

```bash
pip install -r requirements.txt
```

---

## Running the Analysis

1. Clone the repository:

```bash
git clone https://github.com/<your-repo>.git
cd <your-repo>
```

2. Place all required data files in the `Data/` directory.

3. Run the notebooks in sequence:
   - `Codes/modelsALL.ipynb`
   - `Codes/sensitivity.ipynb`

4. Figures will be automatically saved to the `Figures/` folder.

---

## Citation

Rizzetti, N., Fanelli, A., Szewczyk, W.\*, Cescatti, A., Ciscar, J.C., Dubois, G., Ibarreta, D., Figuerola, J., & Massaro, E.  
*Spatial role of land cover on West Nile virus disease in Europe.*  
European Commission, Joint Research Centre (JRC), 2024.

---

## Contact

For scientific questions:  
**Wojciech Szewczyk** — Wojciech.SZEWCZYK@ec.europa.eu  

For code-related queries:  
**Emanuele Massaro**



