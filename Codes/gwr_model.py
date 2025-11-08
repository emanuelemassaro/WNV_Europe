

# --- Imports -----------------------------------------------------------------
import geopandas as gpd
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from libpysal.weights import DistanceBand, lag_spatial
from mgwr.gwr import GWR
from mgwr.sel_bw import Sel_BW
import seaborn as sns
# --- Load shapefile ----------------------------------------------------------
gdf = gpd.read_file('../Data/Shapefiles/20240627_APES FINAL_metric.shp')

# --- Filter & basic cleaning -------------------------------------------------
# remove ZIKV, drop empty geometries
gdf = gdf.loc[gdf['RcrdTyp'] != 'ZIKV'].copy()
gdf = gdf.loc[~gdf.geometry.is_empty & gdf.geometry.notna()].copy()

# columns
id_cols = ["NUTS_NA", "NUTS_CO"]
predictor_cols = [
    "Forest", "Shrub", "Urban", "Crop", "Water", "Other",
    "Mn_P_Wn",   # Mean_P_Winter
    "Mx_WD_W",   # Max_WD_Winter
    "Mx_DD_W",   # Max_DD_Winter
    "Mn_T_Sp",   # Mean_T_Spring
    "Mx_DD_Sp",  # Max_DD_Spring
    "Mn_T_Sm",   # Mean_T_Summer  (included to match your R comment)
    "Mx_WD_A",   # Max_WD_Autumn
    "Mx_DD_A",   # Max_DD_Autumn
    "gdp", "pop_dns"
]
response_col = "incidnc"


# keep only available columns (plus geometry and IDs)
keep = [c for c in (id_cols + predictor_cols + [response_col, "geometry"]) if c in gdf.columns]
gdf = gdf[keep].copy()

# --- Aggregate to NUTS level -------------------------------------------------
# incidence = sum; all predictors = mean (like your R summarise)
agg_dict = {response_col: "sum"}
agg_dict.update({c: "mean" for c in predictor_cols if c in gdf.columns})

fin_eu3 = (
    gdf
    .groupby(id_cols, as_index=False)
    .agg(agg_dict)
)

# restore a representative geometry per NUTS (assumes unique geometry per NUTS)
geom_first = (
    gdf[id_cols + ["geometry"]]
    .drop_duplicates(id_cols)
)
fin_eu3 = fin_eu3.merge(geom_first, on=id_cols, how="left")
fin_eu3 = gpd.GeoDataFrame(fin_eu3, geometry="geometry", crs=gdf.crs)

# drop rows with any NA in variables we need
needed = [response_col] + [c for c in predictor_cols if c in fin_eu3.columns]
fin_eu3 = fin_eu3.dropna(subset=needed).copy()



# --- Project to metric CRS and compute centroids -----------------------------
# Use a metric CRS for distances (UTM zone 33N as in your R code)
fin_eu3 = fin_eu3.to_crs(32633)
centroids = fin_eu3.geometry.centroid  # in meters now
coords = np.c_[centroids.x, centroids.y]



# Replace StandardScaler with manual z-score using ddof=1
num_cols = fin_eu3.select_dtypes(include="number").columns.tolist()

fin_eu3_std = fin_eu3.copy()
for c in num_cols:
    s = fin_eu3[c]
    fin_eu3_std[c] = (s - s.mean()) / s.std(ddof=1)   # R's scale()


y_raw = fin_eu3_std["incidnc"].to_numpy()
X_base = fin_eu3_std[predictor_cols].to_numpy()

kernels = ["gaussian", "bisquare", "exponential"]   # add/remove kernels as you like
distances = list(range(100_000, 200_001, 10_000))   # 100–200 km by 10 km
eps = 1e-9

rows = []
for d in distances:
    # Build binary distance-band weights (no row standardization) → SUM of neighbors
    w = DistanceBand.from_array(coords, threshold=d + eps, binary=True, silence_warnings=True)
    W = w.full()[0]  # dense 0/1 matrix

    # Isolated units diagnostic
    isolated = int((W.sum(axis=1) == 0).sum())

    # Spatial lag of y: Wy = SUM of neighbors' y (binary weights)
    Wy = W.dot(y_raw)
    y = Wy.reshape(-1, 1)
    X = X_base  # predictors unchanged

    for ker in kernels:
        try:
            # adaptive bandwidth by AICc
            bw = Sel_BW(coords, y, X, kernel=ker, fixed=False).search()
            res = GWR(coords, y, X, bw, kernel=ker, fixed=False, constant=True).fit()

            rows.append({
                "kernel": ker,
                "fixed": False,                 # adaptive GWR
                "max_distance_m": d,
                "isolated_units": isolated,
                "bandwidth": bw,                # adaptive: neighbor count
                "AICc": float(res.aicc),
                "mean_localR2": float(res.localR2.mean()),
                "min_localR2":  float(res.localR2.min()),
                "max_localR2":  float(res.localR2.max()),
            })
            print(f"d={d/1000:.0f} km | {ker:10s} | bw={bw:4} | "
                  f"AICc={res.aicc:.1f} | R2 μ/lo/hi = "
                  f"{res.localR2.mean():.3f}/{res.localR2.min():.3f}/{res.localR2.max():.3f} | iso={isolated}")

        except Exception as e:
            rows.append({
                "kernel": ker,
                "fixed": False,
                "max_distance_m": d,
                "isolated_units": isolated,
                "bandwidth": np.nan,
                "AICc": np.nan,
                "mean_localR2": np.nan,
                "min_localR2":  np.nan,
                "max_localR2":  np.nan,
                "error": str(e),
            })
            print(f"d={d/1000:.0f} km | {ker:10s} | ERROR: {e}")

# Save to CSV
out_path = 'sensitivity.csv'
sens_df = pd.DataFrame(rows).sort_values(["max_distance_m", "kernel"]).reset_index(drop=True)
out_path = f"sensitivity_gwr_Wy_binarysum_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"
sens_df.to_csv(out_path, index=False)
print("Saved:", out_path)



