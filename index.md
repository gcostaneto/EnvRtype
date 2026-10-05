# EnvRtype

### Envirotyping for Quantitative Genetics and Plant Breeding

**EnvRtype** is an R package for **enviromics** — the study of the
*envirome*, the set of environmental conditions linked to the biological
performance of living beings (Costa-Neto and Fritsche-Neto, 2021; Crossa
et al., 2021). It collects worldwide daily weather and soil data,
derives agro-meteorological parameters, builds environmental covariable
matrices and relatedness kernels, mines environmental typologies,
delineates soil zones with Gaussian mixture models, and runs a FAO-56
soil water balance for use in enviromics and genotype-by-environment
(GxE) analyses.

This is a fully re-engineered successor to the original [EnvRtype
(allogamous/EnvRtype, 2021)](https://github.com/allogamous/EnvRtype),
extending it from a weather-typing toolbox into an end-to-end
envirotyping, simulation and prediction framework.

> 📄 **Original publication:** Costa-Neto, G., Galli, G., Carvalho, H.
> F., Crossa, J., & Fritsche-Neto, R. (2021). *EnvRtype: a software to
> interplay enviromics and quantitative genomics in agriculture.* **G3
> Genes\|Genomes\|Genetics**, 11(4), jkab040. [Read the paper
> →](https://academic.oup.com/g3journal/article/11/4/jkab040/6129777)

## 📖 **Documentation:** <https://gcostaneto.github.io/EnvRtype/>

## Installation

``` r

if (!require("devtools")) install.packages("devtools")
devtools::install_github("gcostaneto/EnvRtype")
```

``` r

library(EnvRtype)
```

------------------------------------------------------------------------

## What’s new compared to the original EnvRtype (2021)

The original [EnvRtype](https://github.com/allogamous/EnvRtype)
(Costa-Neto et al., 2021, *G3*) focused on collecting NASA POWER weather
data, computing environmental typologies and building environmental
kernels for reaction-norm models. This version keeps that workflow and
adds several new layers:

| Area | Original EnvRtype (2021) | EnvRtype (this version) |
|----|----|----|
| **Weather data** | [`get_weather()`](https://gcostaneto.github.io/EnvRtype/reference/get_weather.md) (NASA POWER, daily) | [`get_weather()`](https://gcostaneto.github.io/EnvRtype/reference/get_weather.md) + hourly ([`get_weather_hourly()`](https://gcostaneto.github.io/EnvRtype/reference/get_weather_hourly.md)) and **resumable / restartable** downloads ([`get_weather_resumable()`](https://gcostaneto.github.io/EnvRtype/reference/get_weather_resumable.md), [`read_progress_log()`](https://gcostaneto.github.io/EnvRtype/reference/read_progress_log.md), [`restart_from_log()`](https://gcostaneto.github.io/EnvRtype/reference/restart_from_log.md)) |
| **Soil data** | — | [`get_soil()`](https://gcostaneto.github.io/EnvRtype/reference/get_soil.md), [`get_soil_resumable()`](https://gcostaneto.github.io/EnvRtype/reference/get_soil_resumable.md), [`soil_classification()`](https://gcostaneto.github.io/EnvRtype/reference/soil_classification.md) (Gaussian-mixture soil zoning) |
| **Other geodata** | — | [`get_elevation()`](https://gcostaneto.github.io/EnvRtype/reference/get_elevation.md), [`get_bioclim()`](https://gcostaneto.github.io/EnvRtype/reference/get_bioclim.md), [`get_spatial()`](https://gcostaneto.github.io/EnvRtype/reference/get_spatial.md), [`get_AEZ()`](https://gcostaneto.github.io/EnvRtype/reference/get_AEZ.md), [`get_climate_scenario()`](https://gcostaneto.github.io/EnvRtype/reference/get_climate_scenario.md) |
| **Processing** | [`processWTH()`](https://gcostaneto.github.io/EnvRtype/reference/processWTH.md), `param_temperature/radiation/atmospheric()`, [`summaryWTH()`](https://gcostaneto.github.io/EnvRtype/reference/summaryWTH.md) | Same, plus a full **FAO-56 water balance** ([`water_balance()`](https://gcostaneto.github.io/EnvRtype/reference/water_balance.md), [`summary_water_balance()`](https://gcostaneto.github.io/EnvRtype/reference/summary_water_balance.md)) and **phenology** ([`env_phenology()`](https://gcostaneto.github.io/EnvRtype/reference/env_phenology.md), [`phenology_templates()`](https://gcostaneto.github.io/EnvRtype/reference/phenology_templates.md), [`planting_window_table()`](https://gcostaneto.github.io/EnvRtype/reference/planting_window_table.md), [`best_planting_date()`](https://gcostaneto.github.io/EnvRtype/reference/best_planting_date.md), [`show_phenology()`](https://gcostaneto.github.io/EnvRtype/reference/phenology_templates.md), [`plot_planting_window()`](https://gcostaneto.github.io/EnvRtype/reference/plot_planting_window.md)) |
| **Characterisation** | [`W_matrix()`](https://gcostaneto.github.io/EnvRtype/reference/W_matrix.md), [`env_typing()`](https://gcostaneto.github.io/EnvRtype/reference/env_typing.md) | Adds [`T_matrix()`](https://gcostaneto.github.io/EnvRtype/reference/T_matrix.md), [`env_indices()`](https://gcostaneto.github.io/EnvRtype/reference/env_indices.md), [`env_expand()`](https://gcostaneto.github.io/EnvRtype/reference/env_expand.md), a full **PCA suite** ([`env_pca()`](https://gcostaneto.github.io/EnvRtype/reference/env_pca.md), [`env_pca_biplot()`](https://gcostaneto.github.io/EnvRtype/reference/env_pca_biplot.md), [`env_pca_scree()`](https://gcostaneto.github.io/EnvRtype/reference/env_pca_scree.md), [`env_loading_curve()`](https://gcostaneto.github.io/EnvRtype/reference/env_loading_curve.md), [`env_pc_associate()`](https://gcostaneto.github.io/EnvRtype/reference/env_pc_associate.md)), correlation tools ([`env_cor()`](https://gcostaneto.github.io/EnvRtype/reference/env_cor.md), [`env_cor_heatmap()`](https://gcostaneto.github.io/EnvRtype/reference/env_cor_heatmap.md)) and coverage/target diagnostics ([`coverage_summary()`](https://gcostaneto.github.io/EnvRtype/reference/coverage_summary.md), [`env_target_importance()`](https://gcostaneto.github.io/EnvRtype/reference/env_target_importance.md)) |
| **Risk & TPE** | — | [`env_risk_profile()`](https://gcostaneto.github.io/EnvRtype/reference/env_risk_profile.md), [`env_copula()`](https://gcostaneto.github.io/EnvRtype/reference/env_copula.md), [`tpe_weights()`](https://gcostaneto.github.io/EnvRtype/reference/tpe_weights.md), [`project_risk()`](https://gcostaneto.github.io/EnvRtype/reference/project_risk.md), [`env_target_importance()`](https://gcostaneto.github.io/EnvRtype/reference/env_target_importance.md) |
| **Kernels** | [`env_kernel()`](https://gcostaneto.github.io/EnvRtype/reference/env_kernel.md), [`get_kernel()`](https://gcostaneto.github.io/EnvRtype/reference/get_kernel.md) | Adds [`decompose_kernels()`](https://gcostaneto.github.io/EnvRtype/reference/decompose_kernels.md), [`undecompose_kernels()`](https://gcostaneto.github.io/EnvRtype/reference/undecompose_kernels.md), [`truncate_gxe_kernel()`](https://gcostaneto.github.io/EnvRtype/reference/truncate_gxe_kernel.md) and a soil kernel path in [`get_kernel()`](https://gcostaneto.github.io/EnvRtype/reference/get_kernel.md) |
| **Modelling** | [`kernel_model()`](https://gcostaneto.github.io/EnvRtype/reference/kernel_model.md) | Adds [`kernel_cv()`](https://gcostaneto.github.io/EnvRtype/reference/kernel_cv.md), [`kernel_model_clustered()`](https://gcostaneto.github.io/EnvRtype/reference/kernel_model_clustered.md), [`kernel_model_mc()`](https://gcostaneto.github.io/EnvRtype/reference/kernel_model_mc.md), [`varcomp_summary()`](https://gcostaneto.github.io/EnvRtype/reference/varcomp_summary.md), environment clustering ([`env_cluster()`](https://gcostaneto.github.io/EnvRtype/reference/env_cluster.md), [`cluster_environments()`](https://gcostaneto.github.io/EnvRtype/reference/cluster_environments.md)) |
| **Untested environments** | — | [`scan_untested_envs()`](https://gcostaneto.github.io/EnvRtype/reference/scan_untested_envs.md), [`grid_scan()`](https://gcostaneto.github.io/EnvRtype/reference/grid_scan.md), [`map_scan()`](https://gcostaneto.github.io/EnvRtype/reference/map_scan.md), [`scan_spatial_table()`](https://gcostaneto.github.io/EnvRtype/reference/scan_spatial_table.md) |
| **Simulation** | — | Ground-truth simulator: [`sim_met()`](https://gcostaneto.github.io/EnvRtype/reference/sim_met.md), [`sim_met_C()`](https://gcostaneto.github.io/EnvRtype/reference/sim_met_C.md), [`sim_W()`](https://gcostaneto.github.io/EnvRtype/reference/sim_W.md), [`sim_W_grid()`](https://gcostaneto.github.io/EnvRtype/reference/sim_W_grid.md) to validate every layer against a known truth |
| **License** | MIT | GPL-3 (CRAN-ready) |

------------------------------------------------------------------------

## Package modules & workflow

The package is organised in four data layers — acquisition, processing,
characterisation and prediction — plus a simulation engine that
generates ground-truth data to validate each layer. Each subsection
lists the functions involved and the primary references used to build
them.

### Overview — the four layers + simulation

*Design follows the enviromics reaction-norm framework of Costa-Neto et
al. (2021, G3) and its envirome-wide extension (Costa-Neto et al., 2023,
G3).*

``` mermaid
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

*Weather and geodata come from NASA POWER via nasapower (Sparks, 2018),
SoilGrids 2.0 (Poggio et al., 2021), WorldClim 2 (Fick & Hijmans, 2017),
bioclimatic predictors (O’Donnell & Ignizio, 2012), CMIP6 / ScenarioMIP
scenarios (Hawkins et al., 2013; O’Neill et al., 2016) and GAEZ v4 (FAO
& IIASA, 2021). Agro-meteorological processing, the FAO-56 soil water
balance and phenology follow Allen et al. (1998), Pereira et al. (2021),
Doorenbos & Kassam (1979), Steduto et al. (2012), Borg & Grimes (1986),
Forsythe et al. (1995), Zadoks et al. (1974) and Counce et al. (2000).*

``` mermaid
flowchart LR
  GW["get_weather()"]
  GWH["get_weather_hourly()"]
  GWR["get_weather_resumable()"]
  RPL["read_progress_log()"]
  RFL["restart_from_log()"]
  GS["get_soil()"]
  GSR["get_soil_resumable()"]
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
  SPH["show_phenology()"]
  PPW["plot_planting_window()"]

  WBAL["water_balance()"]
  SWB["summary_water_balance()"]

  GWR --> RPL --> RFL --> GW
  GW --> PWT
  GWH --> PWT
  PRAD --> PWT
  PTMP --> PWT
  GEL --> PATM
  PATM --> PWT
  PWT --> SWT
  PWT --> EPH
  PHT --> EPH
  EPH --> PWD
  EPH --> SPH
  PWD --> BPD
  PWD --> PPW
  PWT --> WBAL
  GSR --> GS
  GS --> WBAL
  GEL --> WBAL
  WBAL --> SWB

  classDef collect fill:#1565c0,stroke:#0d47a1,color:#fff
  classDef char fill:#2e7d32,stroke:#1b5e20,color:#fff
  classDef wb fill:#0097a7,stroke:#006064,color:#fff
  class GW,GWH,GWR,RPL,RFL,GS,GSR,GEL,GBC,GSP,GAEZ,GCS,PRAD,PATM,PTMP,PWT,SWT collect
  class EPH,PHT,PWD,BPD,SPH,PPW char
  class WBAL,SWB wb
```

### 3 · Characterisation and dissection

*Environmental covariables, typologies and kernels follow Costa-Neto et
al. (2021, G3); interval-resolved indices follow Della Coletta et
al. (2023); breeding-zone typologies and target-population importance
follow Costa-Neto et al. (2023, Agronomy Journal); copula
reparameterisation follows Sklar (1959), Genest & Rivest (1993) and
Salvadori et al. (2007); Gaussian-mixture soil zoning follows Fraley &
Raftery (2002), Scrucca et al. (2016) and Egozcue et al. (2003).*

``` mermaid
flowchart LR
  SWT["summaryWTH()"]
  GS["get_soil()"]
  WM["W_matrix()"]
  ETYP["env_typing()"]
  TM["T_matrix()"]
  EEX["env_expand()"]
  ECOP["env_copula()"]
  ERP["env_risk_profile()"]
  TPW["tpe_weights()"]
  PRK["project_risk()"]
  ETI["env_target_importance()"]
  CVS["coverage_summary()"]
  GCS["get_climate_scenario()"]
  SCL["soil_classification()"]

  EIX["env_indices()"]
  EPC["env_pca()"]
  EPS["env_pca_scree()"]
  EPB["env_pca_biplot()"]
  ELC["env_loading_curve()"]
  ECOR["env_cor()"]
  ECH["env_cor_heatmap()"]
  EPA["env_pc_associate()"]

  SWT --> WM
  WM --> ETYP
  ETYP --> TM
  WM --> EEX
  WM --> ECOP
  WM --> ETI
  ERP --> TPW
  ERP --> PRK
  GCS --> PRK
  TPW --> CVS
  GS --> SCL
  EIX --> WM
  EIX --> EPC
  EIX --> ECOR
  ECOR --> ECH
  EPC --> EPS
  EPC --> EPB
  EPC --> ELC
  EPC --> EPA

  classDef char fill:#2e7d32,stroke:#1b5e20,color:#fff
  classDef dis fill:#6a1b9a,stroke:#4a148c,color:#fff
  classDef collect fill:#1565c0,stroke:#0d47a1,color:#fff
  class WM,ETYP,TM,EEX,ECOP,ERP,TPW,PRK,ETI,CVS,SCL char
  class EIX,EPC,EPS,EPB,ELC,ECOR,ECH,EPA dis
  class SWT,GS,GCS collect
```

### 4 · Prediction and scanning

*Reaction-norm and kernel models follow Costa-Neto et al. (2021, G3),
the nonlinear / Gaussian kernels of Costa-Neto, Fritsche-Neto & Crossa
(2021, Heredity) and the deep kernels of Cuevas et al. (2019); the GxE
modelling engine follows the BGGE approach of Granato et al. (2018);
scanning of untested environments follows the envirome-wide prediction
of Costa-Neto et al. (2023, G3).*

``` mermaid
flowchart TB
  EK["env_kernel()"]
  GK["get_kernel()<br/>K_G · K_E · K_S"]
  DK["decompose_kernels()"]
  UDK["undecompose_kernels()"]
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
  CVS["coverage_summary()"]

  EK --> GK
  GK --> DK
  DK --> UDK
  DK --> TGK
  TGK --> KM
  DK --> KM
  KM --> KCV
  KM --> VCS
  KM --> SUE
  ECL --> CLE
  CLE --> KMC
  KM --> KMM
  SUE --> GSC
  SUE --> MSC
  SUE --> SST
  SUE --> CVS

  classDef pred fill:#c62828,stroke:#b71c1c,color:#fff
  classDef rel fill:#ef6c00,stroke:#e65100,color:#fff
  classDef char fill:#2e7d32,stroke:#1b5e20,color:#fff
  class GK,DK,UDK,TGK,KM,KCV,KMC,KMM,VCS,SUE,GSC,MSC,SST,CVS pred
  class EK rel
  class ECL,CLE char
```

### Simulation — closing the validation loop

*The simulator encodes a known environmental covariance and
reaction-norm structure to benchmark every layer, following the
envirome-wide prediction framework of Costa-Neto et al. (2023, G3) and
the kernel models of Costa-Neto et al. (2021, G3).*

``` mermaid
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

------------------------------------------------------------------------

## Prediction

EnvRtype fits **enviromics-enabled reaction-norm models** in which
phenotypes are regressed on genomic, environmental and GxE relatedness
kernels. The current prediction workflow provides:

- **Kernel construction** —
  [`get_kernel()`](https://gcostaneto.github.io/EnvRtype/reference/get_kernel.md)
  assembles genomic (K_G), environmental (K_E) and soil (K_S) kernels;
  [`env_kernel()`](https://gcostaneto.github.io/EnvRtype/reference/env_kernel.md)
  builds Gaussian / linear environmental kernels;
  [`decompose_kernels()`](https://gcostaneto.github.io/EnvRtype/reference/decompose_kernels.md),
  [`undecompose_kernels()`](https://gcostaneto.github.io/EnvRtype/reference/undecompose_kernels.md)
  and
  [`truncate_gxe_kernel()`](https://gcostaneto.github.io/EnvRtype/reference/truncate_gxe_kernel.md)
  manage the GxE block structure.
- **Model fitting** —
  [`kernel_model()`](https://gcostaneto.github.io/EnvRtype/reference/kernel_model.md)
  fits the multi-kernel Bayesian reaction-norm model;
  [`kernel_model_mc()`](https://gcostaneto.github.io/EnvRtype/reference/kernel_model_mc.md)
  runs multiple chains;
  [`kernel_model_clustered()`](https://gcostaneto.github.io/EnvRtype/reference/kernel_model_clustered.md)
  fits environment-clustered models built with
  [`env_cluster()`](https://gcostaneto.github.io/EnvRtype/reference/env_cluster.md)
  /
  [`cluster_environments()`](https://gcostaneto.github.io/EnvRtype/reference/cluster_environments.md).
- **Validation & variance** —
  [`kernel_cv()`](https://gcostaneto.github.io/EnvRtype/reference/kernel_cv.md)
  runs cross-validation (CV1 / CV2 / CV0 schemes) and
  [`varcomp_summary()`](https://gcostaneto.github.io/EnvRtype/reference/varcomp_summary.md)
  reports variance components and the share of variance explained by
  each kernel.
- **Untested environments** —
  [`scan_untested_envs()`](https://gcostaneto.github.io/EnvRtype/reference/scan_untested_envs.md)
  predicts genotype performance in environments with no phenotypes,
  while
  [`grid_scan()`](https://gcostaneto.github.io/EnvRtype/reference/grid_scan.md),
  [`map_scan()`](https://gcostaneto.github.io/EnvRtype/reference/map_scan.md),
  [`scan_spatial_table()`](https://gcostaneto.github.io/EnvRtype/reference/scan_spatial_table.md)
  and
  [`coverage_summary()`](https://gcostaneto.github.io/EnvRtype/reference/coverage_summary.md)
  summarise the target population of environments (TPE).

These tools implement the reaction-norm and kernel methods of Costa-Neto
et al. (2021, *G3*), the nonlinear / Gaussian kernels of Costa-Neto,
Fritsche-Neto & Crossa (2021, *Heredity*), the deep kernels of Cuevas et
al. (2019, *G3*), the BGGE GxE engine of Granato et al. (2018, *G3*),
and the envirome-wide prediction of untested environments of Costa-Neto
et al. (2023, *G3*). Full citations are listed under **References for
the new functions** below.

------------------------------------------------------------------------

## Acknowledgements

### Original EnvRtype team (2020–2021)

EnvRtype began in 2020 and was first published as
[allogamous/EnvRtype](https://github.com/allogamous/EnvRtype) in
Costa-Neto, G., Galli, G., Carvalho, H. F., Crossa, J., & Fritsche-Neto,
R. (2021). *EnvRtype: a software to interplay enviromics and
quantitative genomics in agriculture.* **G3: Genes\|Genomes\|Genetics**,
11(4), jkab040. We gratefully acknowledge the original authors and the
community of users whose feedback shaped this new version.

### Current maintainers and developers

- **Germano Costa-Neto** — maintainer & lead developer
  \<<germano.cneto@gmail.com>\>
- **Fernanda Pontes** — developer \<<ferpoontes@gmail.com>\>

------------------------------------------------------------------------

## References for the new functions

The methods and data sources underlying the new and extended functions
are documented below, grouped by module.

### Enviromics framework, kernels and genomic prediction

- Costa-Neto, G., Galli, G., Carvalho, H. F., Crossa, J., &
  Fritsche-Neto, R. (2021). EnvRtype: a software to interplay enviromics
  and quantitative genomics in agriculture. *G3:
  Genes\|Genomes\|Genetics*, 11(4), jkab040.
  <https://doi.org/10.1093/g3journal/jkab040>
- Costa-Neto, G., Fritsche-Neto, R., & Crossa, J. (2021). Nonlinear
  kernels, dominance, and envirotyping data increase the accuracy of
  genome-based prediction in multi-environment trials. *Heredity*, 126,
  92–106. <https://doi.org/10.1038/s41437-020-00353-1>
- Costa-Neto, G., Galli, G., Carvalho, H. F., Crossa, J., &
  Fritsche-Neto, R. (2023). Envirome-wide associations enhance
  multi-environment prediction. *G3: Genes\|Genomes\|Genetics*, 13(2),
  jkac313. <https://doi.org/10.1093/g3journal/jkac313>
- Cuevas, J., Montesinos-López, O., Juliana, P., Guzmán, C.,
  Pérez-Rodríguez, P., González-Bucio, J., Burgueño, J.,
  Montesinos-López, A., & Crossa, J. (2019). Deep kernel for genomic and
  near infrared predictions in multi-environment breeding trials. *G3:
  Genes\|Genomes\|Genetics*, 9(9), 2913–2924.
  <https://doi.org/10.1534/g3.119.400493>
- Granato, I., Cuevas, J., Luna-Vázquez, F., Crossa, J.,
  Montesinos-López, O., Burgueño, J., & Fritsche-Neto, R. (2018). BGGE:
  a new package for genomic-enabled prediction incorporating genotype ×
  environment interaction models. *G3: Genes\|Genomes\|Genetics*, 8(9),
  3039–3047. <https://doi.org/10.1534/g3.118.200435>

### Environmental characterisation, typologies and target-population importance

- Costa-Neto, G., da Matta, D., Fernandes, I. K., & Heinemann, A. B.
  (2023). Environmental clusters defining breeding zones for tropical
  irrigated rice in Brazil. *Agronomy Journal*, 115(5).
  <https://doi.org/10.1002/agj2.21481>
- Della Coletta, R., Liese, S. E., Fernandes, S. B., Mikel, M. A.,
  Bohn, M. O., Lipka, A. E., & Hirsch, C. N. (2023). Linking genetic and
  environmental factors through marker effect networks to understand
  trait plasticity. *GENETICS*, 224(4), iyad103.
  <https://doi.org/10.1093/genetics/iyad103>

### Copula-based environmental indices (`env_copula`)

- Sklar, A. (1959). Fonctions de répartition à n dimensions et leurs
  marges. *Publications de l’Institut de Statistique de l’Université de
  Paris*, 8, 229–231.
- Genest, C., & Rivest, L.-P. (1993). Statistical inference procedures
  for bivariate Archimedean copulas. *Journal of the American
  Statistical Association*, 88(423), 1034–1043.
  <https://doi.org/10.1080/01621459.1993.10476372>
- Salvadori, G., De Michele, C., Kottegoda, N. T., & Rosso, R. (2007).
  *Extremes in Nature: An Approach Using Copulas*. Springer.

### Environmental data acquisition (weather, soil, elevation, bioclimate, scenarios, agro-ecological zones)

- Sparks, A. H. (2018). nasapower: NASA-POWER Data from R. *Journal of
  Open Source Software*, 3(30), 1035.
  <https://doi.org/10.21105/joss.01035>
- Forsythe, W. C., Rykiel, E. J., Stahl, R. S., Wu, H., &
  Schoolfield, R. M. (1995). A model comparison for daylength as a
  function of latitude and day of the year. *Ecological Modelling*,
  80(1), 87–95. <https://doi.org/10.1016/0304-3800(94)00034-F>
- Peng, S., Huang, J., Sheehy, J. E., Laza, R. C., Visperas, R. M.,
  Zhong, X., Centeno, G. S., Khush, G. S., & Cassman, K. G. (2004). Rice
  yields decline with higher night temperature from global warming.
  *Proceedings of the National Academy of Sciences*, 101(27), 9971–9975.
  <https://doi.org/10.1073/pnas.0403720101>
- Poggio, L., de Sousa, L. M., Batjes, N. H., Heuvelink, G. B. M.,
  Kempen, B., Ribeiro, E., & Rossiter, D. (2021). SoilGrids 2.0:
  producing soil information for the globe with quantified spatial
  uncertainty. *SOIL*, 7(1), 217–240.
  <https://doi.org/10.5194/soil-7-217-2021>
- Fick, S. E., & Hijmans, R. J. (2017). WorldClim 2: new 1-km spatial
  resolution climate surfaces for global land areas. *International
  Journal of Climatology*, 37(12), 4302–4315.
  <https://doi.org/10.1002/joc.5086>
- O’Donnell, M. S., & Ignizio, D. A. (2012). Bioclimatic predictors for
  supporting ecological applications in the conterminous United States.
  *U.S. Geological Survey Data Series*, 691.
  <https://doi.org/10.3133/ds691>
- Hawkins, E., Osborne, T. M., Ho, C. K., & Challinor, A. J. (2013).
  Calibration and bias correction of climate projections for crop
  modelling: an idealised case study over Europe. *Agricultural and
  Forest Meteorology*, 170, 19–31.
  <https://doi.org/10.1016/j.agrformet.2012.04.007>
- O’Neill, B. C., Tebaldi, C., van Vuuren, D. P., Eyring, V.,
  Friedlingstein, P., Hurtt, G., Knutti, R., Kriegler, E., Lamarque,
  J.-F., Lowe, J., Meehl, G. A., Moss, R., Riahi, K., & Sanderson, B. M.
  (2016). The Scenario Model Intercomparison Project (ScenarioMIP) for
  CMIP6. *Geoscientific Model Development*, 9(9), 3461–3482.
  <https://doi.org/10.5194/gmd-9-3461-2016>
- FAO & IIASA (2021). *Global Agro-Ecological Zones (GAEZ v4)*. FAO,
  Rome. <https://gaez.fao.org>

### Water balance, evapotranspiration and phenology

- Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998). Crop
  evapotranspiration: guidelines for computing crop water requirements.
  *FAO Irrigation and Drainage Paper 56*. FAO, Rome.
- Pereira, L. S., Paredes, P., López-Urrea, R., Hunsaker, D. J., Mota,
  M., & Mohammadi Shad, Z. (2021). Standard single and basal crop
  coefficients for field crops. Updates and advances to the FAO56 crop
  water requirements method. *Agricultural Water Management*,
  243, 106466. <https://doi.org/10.1016/j.agwat.2020.106466>
- Doorenbos, J., & Kassam, A. H. (1979). Yield response to water. *FAO
  Irrigation and Drainage Paper 33*. FAO, Rome.
- Steduto, P., Hsiao, T. C., Fereres, E., & Raes, D. (2012). Crop yield
  response to water. *FAO Irrigation and Drainage Paper 66*. FAO, Rome.
- Borg, H., & Grimes, D. W. (1986). Depth development of roots with
  time: an empirical description. *Transactions of the ASAE*, 29(1),
  194–197. <https://doi.org/10.13031/2013.30125>
- Zadoks, J. C., Chang, T. T., & Konzak, C. F. (1974). A decimal code
  for the growth stages of cereals. *Weed Research*, 14(6), 415–421.
  <https://doi.org/10.1111/j.1365-3180.1974.tb01084.x>
- Counce, P. A., Keisling, T. C., & Mitchell, A. J. (2000). A uniform,
  objective, and adaptive system for expressing rice development. *Crop
  Science*, 40(2), 436–443.
  <https://doi.org/10.2135/cropsci2000.402436x>

### Soil classification and Gaussian-mixture risk zoning

- Fraley, C., & Raftery, A. E. (2002). Model-based clustering,
  discriminant analysis, and density estimation. *Journal of the
  American Statistical Association*, 97(458), 611–631.
  <https://doi.org/10.1198/016214502760047131>
- Scrucca, L., Fop, M., Murphy, T. B., & Raftery, A. E. (2016). mclust
  5: clustering, classification and density estimation using Gaussian
  finite mixture models. *The R Journal*, 8(1), 205–233.
  <https://doi.org/10.32614/RJ-2016-021>
- Egozcue, J. J., Pawlowsky-Glahn, V., Mateu-Figueras, G., &
  Barceló-Vidal, C. (2003). Isometric logratio transformations for
  compositional data analysis. *Mathematical Geology*, 35(3), 279–300.
  <https://doi.org/10.1023/A:1023818214614>
- Zhang, F., Jia, Z., Wu, S., Chen, C., Chen, X., Zheng, C., & Xu, M.
  (2025). Machine learning and Gaussian mixture model for delineating
  soil cadmium risk zones. *Ecosystem Health and Sustainability*.
  <https://doi.org/10.34133/ehs.0402>

------------------------------------------------------------------------

## Citation

If you use **EnvRtype** in your research, please cite the original
publication:

> Costa-Neto, G., Galli, G., Carvalho, H. F., Crossa, J., &
> Fritsche-Neto, R. (2021). *EnvRtype: a software to interplay
> enviromics and quantitative genomics in agriculture.* **G3
> Genes\|Genomes\|Genetics**, 11(4), jkab040.
> <https://doi.org/10.1093/g3journal/jkab040>
>
> 🔗
> [academic.oup.com/g3journal/article/11/4/jkab040/6129777](https://academic.oup.com/g3journal/article/11/4/jkab040/6129777)

BibTeX

``` bibtex
@article{costaneto2021envrtype,
  title   = {EnvRtype: a software to interplay enviromics and quantitative genomics in agriculture},
  author  = {Costa-Neto, Germano and Galli, Giovanni and Carvalho, Humberto Fanelli
             and Crossa, Jos{\'e} and Fritsche-Neto, Roberto},
  journal = {G3 Genes|Genomes|Genetics},
  volume  = {11},
  number  = {4},
  pages   = {jkab040},
  year    = {2021},
  doi     = {10.1093/g3journal/jkab040},
  url     = {https://academic.oup.com/g3journal/article/11/4/jkab040/6129777}
}
```

From within R:

``` r

citation("EnvRtype")
```

Please also kindly let the maintainer know when you are using EnvRtype
(contact **<germano.cneto@gmail.com>**).

------------------------------------------------------------------------

## License

Released under the **GPL-3** license — free and open source. You may
use, study, modify and redistribute the code, provided derivative works
remain under the same terms.
