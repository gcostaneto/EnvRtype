# Estimate Thermal Parameters

Computes growing degree days, the temperature-limited radiation use
efficiency factor, and daily temperature range.

## Usage

``` r
param_temperature(
  env.data,
  Tmax = NULL,
  Tmin = NULL,
  Tbase1 = 9,
  Tbase2 = 45,
  Topt1 = 26,
  Topt2 = 32,
  method = c("standard", "legacy"),
  merge = FALSE,
  verbose = TRUE
)
```

## Arguments

- env.data:

  data.frame. A `get_weather` output or equivalent.

- Tmax, Tmin:

  character. Column names for maximum and minimum air temperature
  (Celsius).

- Tbase1:

  numeric. Base temperature below which development stops.

- Tbase2:

  numeric. Upper temperature at which development stops.

- Topt1, Topt2:

  numeric. Lower and upper bounds of the optimal range.

- method:

  character. GDD formulation: `"standard"` (default) or `"legacy"`. See
  Details.

- merge:

  boolean. Bind results onto `env.data`.

- verbose:

  boolean. Print a summary of what was computed.

## Value

A data.frame with `GDD` (C/day), `FRUE` (0-1) and `T2M_RANGE` (C).

## Details

**GDD methods.**

- `"standard"`:

  \\GDD = \max(\frac{\min(T\_{max}, T\_{base2}) + T\_{min}}{2} -
  T\_{base1},\\ 0)\\. Tmax is capped at the upper threshold; the result
  is floored at zero. This is the conventional agronomic definition.

- `"legacy"`:

  Clamps *both* Tmax and Tmin into \\\[T\_{base1}, T\_{base2}\]\\ before
  averaging. This is the behaviour of earlier EnvRtype versions and
  corresponds to the "Method 2" variant. It credits heat units on days
  colder than \\T\_{base1}\\ – e.g. with Tmax 12, Tmin 5, Tbase1 9 it
  returns 1.5 where the standard method returns 0.

The default changed from legacy to standard; results will differ on cold
days. Pass `method = "legacy"` to reproduce earlier output.

`FRUE` rises linearly from 0 at `Tbase1` to 1 at `Topt1`, holds at 1
across the optimum, and falls to 0 at `Tbase2`.

## See also

`param_radiation`, `param_atmospheric`, `processWTH`

## Author

Germano Costa Neto

## Examples

``` r
if (FALSE) { # \dontrun{
env.data <- get_weather(lat = -13.05, lon = -56.05)
param_temperature(env.data)
param_temperature(env.data, method = "legacy")  # pre-refactor behaviour
} # }
```
