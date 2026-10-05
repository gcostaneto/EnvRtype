# Daily Soil Water Balance and Crop Water-Stress Indices

Runs a FAO-56 single crop coefficient soil water balance by combining
daily weather from `get_weather` with soil hydraulic properties from
`get_soil`, returning daily root-zone depletion, a water-stress
coefficient, actual evapotranspiration, and a categorical stress level.

The stress coefficient `Ks` is the main output: bounded in \[0, 1\],
dimensionless, and directly proportional to the reduction in
transpiration.

## Usage

``` r
water_balance(
  env.data,
  soil.data,
  env.id = "env",
  days.id = "daysFromStart",
  PREC = "PRECTOT",
  ETo = NULL,
  irrigation = 0,
  root.depth = 1000,
  root.init = 300,
  dap.root.max = 60,
  kc.stages = c(1, 20, 45, 90, 120),
  kc.values = c(0.4, 0.75, 1.15, 1.15, 0.6),
  p = 0.55,
  p.adjust = TRUE,
  initial.depletion = NULL,
  fc.var = "wv0033",
  pwp.var = "wv1500",
  convert.units = TRUE,
  verbose = TRUE
)
```

## Arguments

- env.data:

  data.frame. Daily weather, typically from `get_weather` or
  `processWTH`. Must contain the environment id, a day counter,
  precipitation, and a reference evapotranspiration column.

- soil.data:

  data.frame. Soil properties from `get_soil(wide = TRUE)`. Must include
  `wv0033` (field capacity) and `wv1500` (wilting point) for the layers
  spanning the root zone.

- env.id:

  character. Environment-id column present in both inputs. Default
  `"env"`.

- days.id:

  character. Column giving days from planting. Default
  `"daysFromStart"`.

- PREC:

  character. Precipitation column (mm/day). Default `"PRECTOT"`.

- ETo:

  character. Reference evapotranspiration column (mm/day). If `NULL`
  (default), the function looks for `"ETP"` (added by
  `param_atmospheric`), then `"EVPTRNS"`.

- irrigation:

  numeric or character. Irrigation in mm/day: a single value applied
  daily, a vector as long as `env.data`, or the name of a column.
  Default 0.

- root.depth:

  numeric. Maximum rooting depth in mm. Default 1000.

- root.init:

  numeric. Rooting depth at emergence in mm. Default 300.

- dap.root.max:

  numeric. Days after planting at which maximum rooting depth is
  reached. Default 60.

- kc.stages, kc.values:

  numeric. Days after planting and the corresponding crop coefficients,
  linearly interpolated between. Defaults describe a generic 120-day
  maize crop (FAO-56 Table 12). **Calibrate these.**

- p:

  numeric. Soil-water depletion fraction for no stress (FAO-56 Table
  22). Default 0.55 for maize.

- p.adjust:

  logical. Adjust `p` daily for evaporative demand using FAO-56 Eq. 84.
  Default `TRUE`.

- initial.depletion:

  numeric. Root-zone depletion at planting in mm, or `NULL` (default) to
  start at field capacity.

- fc.var, pwp.var:

  character. SoilGrids property names for field capacity and wilting
  point. Defaults `"wv0033"` and `"wv1500"`.

- convert.units:

  logical. Apply ISRIC conversion factors to the soil data. Set `FALSE`
  only if already converted. Default `TRUE`.

- verbose:

  logical. Print progress messages. Default `TRUE`.

## Value

A data.frame with one row per environment-day, containing the input
columns plus:

- `TAW`:

  total available water in the root zone (mm)

- `RAW`:

  readily available water, \\p \times TAW\\ (mm)

- `Dr`:

  root-zone depletion at end of day (mm)

- `ASW`:

  available soil water remaining, \\TAW - D_r\\ (mm)

- `depletion_frac`:

  \\D_r / TAW\\; 0 at field capacity, 1 at wilting point

- `Ks`:

  water-stress coefficient; 1 unstressed, 0 fully stressed

- `ETc`:

  crop evapotranspiration, standard conditions (mm/day)

- `ETa`:

  actual evapotranspiration, \\K_s \times ET_c\\

- `deficit`:

  \\ET_c - ET_a\\ (mm/day)

- `drainage`:

  deep percolation below the root zone (mm/day)

- `irrigation`:

  irrigation applied (mm/day)

- `Kc`:

  interpolated crop coefficient

- `root_depth`:

  rooting depth that day (mm)

- `stress_level`:

  ordered factor: none, mild, moderate, severe

## Details

**Method.** FAO-56 single crop coefficient approach (Allen et al. 1998,
Ch. 8). Root-zone depletion is stepped forward daily: \$\$D\_{r,i} =
D\_{r,i-1} - P_i - I_i + ET\_{a,i} + DP_i\$\$ bounded to \\\[0, TAW\]\\.
Water in excess of field capacity becomes deep percolation; depletion
cannot exceed total available water.

The stress coefficient follows FAO-56 Eq. 84: \$\$K_s = 1 \quad
\mathrm{if}\\ D_r \le RAW\$\$ \$\$K_s = \frac{TAW - D_r}{TAW - RAW}
\quad \mathrm{if}\\ D_r \> RAW\$\$ so transpiration proceeds at the
potential rate until readily available water is exhausted, then declines
linearly to zero at wilting point.

**Causal ordering.** `Ks` on day \\i\\ is computed from *yesterday's*
depletion against *today's* capacity. Using the same day's closing
depletion would let today's rainfall relieve stress that already reduced
today's transpiration – a subtle look-ahead that systematically
understates stress.

**Units.** SoilGrids serves integers in mapped units: `wv0033 = 314`
means 0.314 cm3/cm3. Conversion is automatic; see `convert.units`.

**Total available water** is integrated across the SoilGrids layers
intersecting the root zone, each weighted by the depth it contributes,
not taken from a single layer. Because rooting depth grows over the
season, TAW is recomputed daily (cached by unique depth for speed).

**What this does not model.** No dual crop coefficient, so soil
evaporation is not separated from transpiration. No capillary rise, no
lateral flow, and **no surface runoff** – all rainfall is assumed to
infiltrate, which overestimates stored water and understates stress on
steep, crusted or low-infiltration soils. Subtract runoff from the
precipitation column beforehand if that matters. Irrigation is applied
as given; the balance does not schedule it.

**Interpreting the output.** Regress `Ks` against yield: it is bounded,
dimensionless, and scales transpiration directly. `depletion_frac` is
easier to explain to agronomists. `stress_level` is a convenience for
grouping and tabulation and should not be used as a modelling covariate
– it discards magnitude.

## References

Allen, R.G., Pereira, L.S., Raes, D. & Smith, M. (1998) Crop
evapotranspiration: guidelines for computing crop water requirements.
*FAO Irrigation and Drainage Paper 56*. FAO, Rome.

Pereira, L.S. et al. (2021) Standard single and basal crop coefficients
for field crops. *Agricultural Water Management* 243, 106466.
[doi:10.1016/j.agwat.2020.106466](https://doi.org/10.1016/j.agwat.2020.106466)

Doorenbos, J. & Kassam, A.H. (1979) Yield response to water. *FAO
Irrigation and Drainage Paper 33*. FAO, Rome.

Steduto, P., Hsiao, T.C., Fereres, E. & Raes, D. (2012) Crop yield
response to water. *FAO Irrigation and Drainage Paper 66*. FAO, Rome.

Poggio, L. et al. (2021) SoilGrids 2.0. *SOIL* 7, 217-240.
[doi:10.5194/soil-7-217-2021](https://doi.org/10.5194/soil-7-217-2021)

Turek, M.E. et al. (2023) Global mapping of volumetric water retention
at 100, 330 and 15000 cm suction using the WoSIS database.
*International Soil and Water Conservation Research* 11(2), 225-239.
[doi:10.1016/j.iswcr.2022.08.001](https://doi.org/10.1016/j.iswcr.2022.08.001)

Borg, H. & Grimes, D.W. (1986) Depth development of roots with time.
*Transactions of the ASAE* 29(1), 194-197.
[doi:10.13031/2013.30125](https://doi.org/10.13031/2013.30125)

## See also

`summary_water_balance` to condense the daily output, `get_weather` and
`get_soil` for the inputs, `processWTH` to derive ETo,
`soil_classification` for cross-site soil zoning.

## Examples

``` r
if (FALSE) { # \dontrun{
## ---------------------------------------------------------------
## 1. Minimal use
## ---------------------------------------------------------------
wth  <- get_weather(env.id = "NM", lat = -13.05, lon = -56.05,
                    start.day = "2015-02-15", end.day = "2015-06-15")
wth  <- processWTH(wth)                 # adds ETP
soil <- get_soil(env.id = "NM", lat = -13.05, lon = -56.05,
                 variables.names = c("wv0033", "wv1500", "clay"))

wb <- water_balance(wth, soil)
head(wb[, c("env", "daysFromStart", "TAW", "Dr", "Ks",
            "deficit", "stress_level")])


## ---------------------------------------------------------------
## 2. Several environments at once
## ---------------------------------------------------------------
sites <- data.frame(
  env = c("SOR", "LON", "BAR"),
  lat = c(-12.5453, -23.3045, -12.1530),
  lon = c(-55.7113, -51.1696, -44.9900))

wth  <- processWTH(get_weather(env.id = sites$env, lat = sites$lat,
                               lon = sites$lon,
                               start.day = "2023-10-15",
                               end.day   = "2024-02-12"))
soil <- get_soil(env.id = sites$env, lat = sites$lat, lon = sites$lon,
                 variables.names = c("wv0033", "wv1500"))

wb <- water_balance(wth, soil, root.depth = 1000)
summary_water_balance(wb)


## ---------------------------------------------------------------
## 3. Irrigation: constant, per-day, or from a column
## ---------------------------------------------------------------
water_balance(wth, soil, irrigation = 5)              # 5 mm every day
water_balance(wth, soil, irrigation = rep(0, nrow(wth)))

wth$irr <- ifelse(wth$daysFromStart %in% 40:70, 8, 0)  # top-up at flowering
wb_irr  <- water_balance(wth, soil, irrigation = "irr")

## How much stress did irrigation remove?
c(rainfed   = sum(water_balance(wth, soil, verbose = FALSE)$deficit),
  irrigated = sum(wb_irr$deficit))


## ---------------------------------------------------------------
## 4. A different crop: soybean
## ---------------------------------------------------------------
## Kc and p come from FAO-56 Tables 12 and 22. Always check them
## against your own crop, cycle length and region.
wb_soy <- water_balance(
  wth, soil,
  kc.stages = c(1, 25, 50, 95, 120),
  kc.values = c(0.40, 0.80, 1.15, 1.15, 0.50),
  p         = 0.50,          # soybean, FAO-56 Table 22
  root.depth = 800)


## ---------------------------------------------------------------
## 5. Restricted root zone (compaction, shallow profile)
## ---------------------------------------------------------------
shallow <- water_balance(wth, soil, root.depth = 400, verbose = FALSE)
deep    <- water_balance(wth, soil, root.depth = 1200, verbose = FALSE)
c(shallow_TAW = max(shallow$TAW), deep_TAW = max(deep$TAW))
## A shallower profile stores less water and stresses sooner.


## ---------------------------------------------------------------
## 6. Dry planting: start part-depleted rather than at field capacity
## ---------------------------------------------------------------
wb_dry <- water_balance(wth, soil, initial.depletion = 60)


## ---------------------------------------------------------------
## 7. Season and growth-stage summaries
## ---------------------------------------------------------------
summary_water_balance(wb)

summary_water_balance(wb, by.interval = TRUE,
                      time.window  = c(0, 14, 35, 60, 90, 120),
                      names.window = c("P-E", "E-V1", "V1-V4",
                                       "V4-VT", "VT-GF", "GF-PM"))
## Stress at VT-GF (flowering to grain fill) costs far more yield
## than the same stress during vegetative growth.


## ---------------------------------------------------------------
## 8. Using the indices as envirotype covariates
## ---------------------------------------------------------------
W <- summary_water_balance(wb, by.interval = TRUE,
                           time.window = c(0, 30, 60, 90, 120))
wide <- reshape(W[, c("env", "interval", "mean_Ks")],
                idvar = "env", timevar = "interval", direction = "wide")
head(wide)   # one row per environment, ready to join to trial data
} # }
```
