# EnvRtype <img src="https://img.shields.io/badge/version-0.1.0-blue" align="right"/>

### Envirotyping for Quantitative Genetics and Plant Breeding

**EnvRtype** is an R package for **enviromics** — the study of the *envirome*, the set of
environmental conditions linked to the biological performance of living beings. It collects
worldwide daily weather and soil data, derives agro-meteorological parameters, builds
environmental covariable matrices and relatedness kernels, mines environmental typologies,
delineates soil zones with Gaussian mixture models, and runs a FAO-56 soil water balance for
use in enviromics and genotype-by-environment (GxE) analyses.

This is a fully re-engineered successor to the original
[EnvRtype (allogamous/EnvRtype, 2021)](https://github.com/allogamous/EnvRtype), extending it
from a weather-typing toolbox into an end-to-end envirotyping, simulation and prediction
framework.

---

## Installation

```r
if (!require("devtools")) install.packages("devtools")
devtools::install_github("gcostaneto/EnvRtype")
```

```r
library(EnvRtype)
```

---

## What's new compared to the original EnvRtype (2021)

The original [EnvRtype](https://github.com/allogamous/EnvRtype) (Costa-Neto et al., 2021,
*G3*) focused on collecting NASA POWER weather data, computing environmental typologies
and building environmental kernels for reaction-norm models. This version keeps that workflow and
adds several new layers:

| Area | Original EnvRtype (2021) | EnvRtype (this version) |
|------|--------------------------|-------------------------|
| **Weather data** | `get_weather()` (NASA POWER, daily) | `get_weather()` + hourly (`get_weather_hourly()`) and **resumable / restartable** downloads (`get_weather_resumable()`, `read_progress_log()`, `restart_from_log()`) |
| **Soil data** | — | `get_soil()`, `get_soil_resumable()`, `soil_classification()` (Gaussian-mixture soil zoning) |
| **Other geodata** | — | `get_elevation()`, `get_bioclim()`, `get_spatial()`, `get_AEZ()`, `get_climate_scenario()` |
| **Processing** | `processWTH()`, `param_temperature/radiation/atmospheric()`, `summaryWTH()` | Same, plus a full **FAO-56 water balance** (`water_balance()`, `summary_water_balance()`) and **phenology** (`env_phenology()`, `phenology_templates()`, `planting_window_table()`, `best_planting_date()`) |
| **Characterisation** | `W_matrix()`, `env_typing()` | Adds `T_matrix()`, `env_indices()`, a full **PCA suite** (`env_pca()`, `env_pca_biplot()`, `env_pca_scree()`, `env_loading_curve()`, `env_pc_associate()`), correlation tools (`env_cor()`, `env_cor_heatmap()`) and coverage/target diagnostics |
| **Risk & TPE** | — | `env_risk_profile()`, `env_copula()`, `tpe_weights()`, `project_risk()`, `env_target_importance()` |
| **Kernels** | `env_kernel()`, `get_kernel()` | Adds `decompose_kernels()`, `undecompose_kernels()`, `truncate_gxe_kernel()` and a soil kernel path in `get_kernel()` |
| **Modelling** | `kernel_model()` | Adds `kernel_cv()`, `kernel_model_clustered()`, `kernel_model_mc()`, `varcomp_summary()`, environment clustering (`env_cluster()`, `cluster_environments()`) |
| **Untested environments** | — | `scan_untested_envs()`, `grid_scan()`, `map_scan()`, `scan_spatial_table()` |
| **Simulation** | — | Ground-truth simulator: `sim_met()`, `sim_met_C()`, `sim_W()`, `sim_W_grid()` to validate every layer against a known truth |
| **License** | MIT | GPL-3 (CRAN-ready) |

---

## Package modules & workflow

The package is organised in five layers, plus a simulation engine that generates ground-truth
data to validate each layer.

### Overview — the five layers

```mermaid
flowchart TB
  L1["<b>1 · ACQUIRE</b><br/>get_weather · get_soil · get_elevation<br/>get_bioclim · get_spatial · get_AEZ<br/>get_climate_scenario · + resumable/hourly variants"]
  L2["<b>2 · PROCESS</b><br/>processWTH · summaryWTH · param_*<br/>water_balance · env_phenology"]
  L3["<b>3 · CHARACTERISE</b><br/>W_matrix · T_matrix · env_indices · env_pca<br/>soil_classification · env_typing · env_risk_profile"]
  L4["<b>4 · MODEL &amp; SCAN</b><br/>env_kernel · get_kernel · kernel_model<br/>kernel_cv · scan_untested_envs · map_scan"]
  SIM["<b>SIMULATE</b><br/>sim_met() → C_env → sim_W()<br/><i>ground truth for every layer</i>"]

  L1 --> L2 --> L3 --> L4
  SIM -.->|"validates"| L3
  SIM -.->|"validates"| L4

  classDef l1 fill:#1565c0,stroke:#0d47a1,color:#fff
  classDef l2 fill:#1976d2,stroke:#0d47a1,color:#fff
  classDef l3 fill:#2e7d32,stroke:#1b5e20,color:#fff
  classDef l4 fill:#c62828,stroke:#b71c1c,color:#fff
  classDef sim fill:#6a1b9a,stroke:#4a148c,color:#fff
  class L1 l1
  class L2 l2
  class L3 l3
  class L4 l4
  class SIM sim
```

### 1–2 · Data acquisition and processing

```mermaid
flowchart LR
  GW["get_weather()"]
  GWH["get_weather_hourly()"]
  GS["get_soil()"]
  GEL["get_elevation()"]
  GBC["get_bioclim()"]
  GSP["get_spatial()"]
  GAEZ["get_AEZ()"]
  GCS["get_climate_scenario()"]

  PRAD["param_radiation()"]
  PATM["param_atmospheric()"]
  PTMP["param_temperature()"]
  PWT["processWTH()"]
  SWT["summaryWTH()"]

  EPH["env_phenology()"]
  PHT["phenology_templates()"]
  PWD["planting_window_table()"]
  BPD["best_planting_date()"]

  WBAL["water_balance()"]
  SWB["summary_water_balance()"]

  GW --> PWT
  PRAD --> PWT
  PTMP --> PWT
  GEL --> PATM
  PATM --> PWT
  PWT --> SWT
  PWT --> EPH
  PHT --> EPH
  EPH --> PWD
  PWD --> BPD
  PWT --> WBAL
  GS --> WBAL
  GEL --> WBAL
  WBAL --> SWB

  classDef collect fill:#1565c0,stroke:#0d47a1,color:#fff
  classDef char fill:#2e7d32,stroke:#1b5e20,color:#fff
  classDef wb fill:#0097a7,stroke:#006064,color:#fff
  class GW,GWH,GS,GEL,GBC,GSP,GAEZ,GCS,PRAD,PATM,PTMP,PWT,SWT collect
  class EPH,PHT,PWD,BPD char
  class WBAL,SWB wb
```

### 3 · Characterisation and dissection

```mermaid
flowchart LR
  SWT["summaryWTH()"]
  WM["W_matrix()"]
  ETYP["env_typing()"]
  TM["T_matrix()"]
  ECOP["env_copula()"]
  ERP["env_risk_profile()"]
  TPW["tpe_weights()"]
  PRK["project_risk()"]
  GCS["get_climate_scenario()"]
  SCL["soil_classification()"]

  EIX["env_indices()"]
  EPC["env_pca()"]
  EPS["env_pca_scree()"]
  EPB["env_pca_biplot()"]
  ELC["env_loading_curve()"]
  ECH["env_cor_heatmap()"]
  EPA["env_pc_associate()"]

  SWT --> WM
  ETYP --> TM
  WM --> ECOP
  ERP --> TPW
  ERP --> PRK
  GCS --> PRK
  EIX --> WM
  EIX --> EPC
  EIX --> ECH
  EPC --> EPS
  EPC --> EPB
  EPC --> ELC
  EPC --> EPA

  classDef char fill:#2e7d32,stroke:#1b5e20,color:#fff
  classDef dis fill:#6a1b9a,stroke:#4a148c,color:#fff
  classDef collect fill:#1565c0,stroke:#0d47a1,color:#fff
  class WM,ETYP,TM,ECOP,ERP,TPW,PRK,SCL char
  class EIX,EPC,EPS,EPB,ELC,ECH,EPA dis
  class SWT,GCS collect
```

### 4 · Prediction and scanning

```mermaid
flowchart TB
  EK["env_kernel()"]
  GK["get_kernel()<br/>K_G · K_E · K_S"]
  DK["decompose_kernels()"]
  TGK["truncate_gxe_kernel()"]
  KM["kernel_model()"]
  KCV["kernel_cv()"]
  KMC["kernel_model_clustered()"]
  KMM["kernel_model_mc()"]
  ECL["env_cluster()"]
  CLE["cluster_environments()"]
  VCS["varcomp_summary()"]
  SUE["scan_untested_envs()"]
  GSC["grid_scan()"]
  MSC["map_scan()"]
  SST["scan_spatial_table()"]

  EK --> GK
  GK --> DK
  DK --> TGK
  TGK --> KM
  DK --> KM
  KM --> KCV
  KM --> VCS
  KM --> SUE
  ECL --> KMC
  CLE --> KMC
  KM --> KMM
  SUE --> GSC
  SUE --> MSC
  SUE --> SST

  classDef pred fill:#c62828,stroke:#b71c1c,color:#fff
  classDef rel fill:#ef6c00,stroke:#e65100,color:#fff
  classDef char fill:#2e7d32,stroke:#1b5e20,color:#fff
  class GK,DK,TGK,KM,KCV,KMC,KMM,VCS,SUE,GSC,MSC,SST pred
  class EK rel
  class ECL,CLE char
```

### Simulation — closing the validation loop

```mermaid
flowchart LR
  K["genomic kinship K<br/><i>n × n lines</i>"]
  SMC["sim_met_C()<br/><i>q × q envs</i>"]
  SM["sim_met()"]
  CE["$C_env<br/><i>q × q correlation<br/>among environments</i>"]
  PH["$data<br/>env · gid · rep · value"]
  SW["sim_W()"]
  SWG["sim_W_grid()"]
  W["simulated W<br/><i>known % explained</i>"]
  GK["get_kernel()"]
  KM["kernel_model()"]
  VCS["varcomp_summary()"]
  CHK{{"compare to<br/>known truth"}}

  K --> SM
  SMC --> SM
  SM --> CE
  SM --> PH
  CE --> SW
  SW --> W
  SWG --> W
  W --> GK
  PH --> KM
  GK --> KM
  KM --> VCS
  VCS --> CHK
  SM -.->|"$truth"| CHK

  classDef sim fill:#6a1b9a,stroke:#4a148c,color:#fff
  classDef pred fill:#c62828,stroke:#b71c1c,color:#fff
  classDef data fill:#455a64,stroke:#263238,color:#fff
  classDef chk fill:#f9a825,stroke:#f57f17,color:#000
  class SM,SW,SWG,SMC sim
  class GK,KM,VCS pred
  class K,CE,PH,W data
  class CHK chk
```

---

## Acknowledgements

### Original EnvRtype team (2020–2021)

EnvRtype began in 2020 and was first published as
[allogamous/EnvRtype](https://github.com/allogamous/EnvRtype) in Costa-Neto, G., Galli, G.,
Carvalho, H. F., Crossa, J., & Fritsche-Neto, R. (2021).
*EnvRtype: a software to interplay enviromics and quantitative genomics in agriculture.*
**G3: Genes|Genomes|Genetics**, 11(4), jkab040. We gratefully acknowledge the original authors
and the community of users whose feedback shaped this new version.

### Current maintainers and developers

- **Germano Costa-Neto** — maintainer & lead developer &lt;germano.cneto@gmail.com&gt;
- **Fernanda Pontes** — developer &lt;ferpoontes@gmail.com&gt;

---

## Citation

```r
citation("EnvRtype")
```

Please cite the package when using it in publications, and kindly let the maintainer know when
you are using EnvRtype (contact **germano.cneto@gmail.com**).

---

## License

Released under the **GPL-3** license — free and open source. You may use, study, modify and
redistribute the code, provided derivative works remain under the same terms.
