# TROLL WTD / Capillary Rise Sensitivity Tests

This folder contains code versions and model outputs used to diagnose the implementation of vertical water movement and water table depth (WTD) effects in TROLL.

The main goal of these tests is to understand whether the simulated soil water dynamics are controlled by physically plausible hydraulic disconnection from the water table, or by numerical artifacts introduced by safeguards such as `theta_w_min`, `Ks_floor`, interface conductivity averaging, vertical discretization, or water-content clamps.

---

## 1. Background

The original TROLL soil-water scheme was based mostly on a bucket-type vertical redistribution logic. The current implementation adds a Darcy-like vertical flux between adjacent soil layers, allowing both upward and downward water movement depending on the hydraulic potential gradient.

The vertical flux is computed from:

- soil water content (`SWC`),
- relative water content / effective saturation (`theta_w`),
- matric potential / soil water potential (`SWP`),
- hydraulic conductivity (`Ks`),
- interface conductivity between adjacent layers,
- vertical hydraulic gradient,
- donor and receiver capacity constraints.

The implementation is intended to represent capillary rise from the water table and vertical redistribution through the soil column, while preventing nonphysical water creation or destruction.

---

## 2. Main diagnostic question

The central question is:

> Is the simulated hydraulic disconnection between the water table and the upper soil layers physically plausible, or is it caused by numerical choices such as `Ks_floor`, `theta_w_min`, interface conductivity averaging, timestep length, or layer thickness?

This is especially important for the deep water table + sandy soil case, where the model often shows strong exchange at the deepest interface but little propagation to upper layers.

---

## 3. Current code assumptions for sensitivity tests

Unless otherwise stated, the tests in this folder use:

```cpp
theta_w_min = 1e-3;
```

This value is used to avoid extremely small or negative `theta_w` values that can generate unstable or undefined hydraulic properties.

The interface hydraulic conductivity is computed using the geometric mean:

```cpp
double k1d = static_cast<double>(Ks_cap[l][d]);
double k2d = static_cast<double>(Ks_cap[l+1][d]);

if (std::isfinite(k1d) && std::isfinite(k2d) && k1d > 0.0 && k2d > 0.0) {
    Ks_cap_interface[l][d] = static_cast<float>(std::sqrt(k1d * k2d));
} else {
    Ks_cap_interface[l][d] = 0.0f;
}
```

Older versions used a threshold such as:

```cpp
if (sum_k > 2e-11f)
```

This threshold should not be used in the current sensitivity tests because it can artificially set small but positive interface conductivities to zero, especially when testing low `Ks_floor` values such as `1e-16`, `1e-18`, or `1e-20`.

If the code still uses the variable name `Ks_cap_harmonic`, note that this name may be outdated. In the current geometric-mean tests, this variable represents interface conductivity, not harmonic conductivity. A clearer name would be:

```cpp
Ks_cap_interface
```

---

## 4. Parameters being tested

### 4.1 `theta_w_min`

Tested values:

```cpp
theta_w_min = 1e-3;
theta_w_min = 1e-6;
```

Preliminary interpretation:

- Changing `theta_w_min` from `1e-3` to `1e-6` did not strongly change the long-term behavior.
- Therefore, `theta_w_min` does not appear to be the main control of the current behavior.
- For the current sensitivity tests, `theta_w_min = 1e-3` is kept fixed.

### 4.2 `Ks_floor`

Tested and proposed values:

```cpp
Ks_floor = 1e-10;
Ks_floor = 1e-12;
Ks_floor = 1e-14;
Ks_floor = 1e-16;
Ks_floor = 1e-18;
Ks_floor = 1e-20;
```

Preliminary interpretation:

- `Ks_floor = 1e-10` appears to keep dry layers too conductive. This can create large vertical fluxes when `redprec` generates strong matric-potential gradients.
- `Ks_floor = 1e-20` appears more stable in terms of SWP, but may strongly limit propagation above the deepest interface.
- `Ks_floor = 1e-16` with geometric mean still appears to produce strong hydraulic disconnection above interface `4_5`.
- Intermediate values such as `1e-12` and `1e-14` should be tested as possible compromises.

---

## 5. Experimental design

### 5.1 First sensitivity block: deepWT sandy

Start with the most sensitive case:

```text
vegetation = veg
water_table = deepWT
texture = sandy
climate = regclim and redprec
```

Run the following configurations:

| Test name | theta_w_min | Ks_floor | Interface K |
|---|---:|---:|---|
| theta1e-3_ks1e-10_geom | 1e-3 | 1e-10 | geometric |
| theta1e-3_ks1e-12_geom | 1e-3 | 1e-12 | geometric |
| theta1e-3_ks1e-14_geom | 1e-3 | 1e-14 | geometric |
| theta1e-3_ks1e-16_geom | 1e-3 | 1e-16 | geometric |
| theta1e-3_ks1e-18_geom | 1e-3 | 1e-18 | geometric |
| theta1e-3_ks1e-20_geom | 1e-3 | 1e-20 | geometric |

The priority tests are:

```text
Ks_floor = 1e-12
Ks_floor = 1e-14
Ks_floor = 1e-16
```

because previous results suggest that `1e-10` may be too conductive and `1e-20` may be too restrictive.

### 5.2 Second sensitivity block: expand to WTD × texture

After selecting two or three candidate `Ks_floor` values, expand to:

```text
deepWT_sandy
deepWT_clayey
shallowWT_sandy
shallowWT_clayey
```

for both:

```text
regclim
redprec
```

Expected qualitative behavior:

- `deepWT_sandy`: strongest hydraulic disconnection from the water table.
- `shallowWT_sandy`: stronger water-table buffering than deepWT.
- `clayey`: higher water retention, lower conductivity, potentially different capillary behavior.
- `redprec`: drier upper layers than `regclim`, especially in deepWT.

### 5.3 Third sensitivity block: geometric vs harmonic interface K

After selecting a candidate `Ks_floor`, compare geometric and harmonic interface conductivity.

Geometric mean:

```cpp
Ks_interface = sqrt(k1 * k2);
```

Harmonic mean:

```cpp
Ks_interface = (2.0 * k1 * k2) / (k1 + k2);
```

This test should be used diagnostically. If the two methods produce completely different regimes, the system is highly sensitive to interface conductivity and the vertical discretization may need further attention.

---

## 6. Key output files

The most important diagnostic files are:

```text
(null)_0_vertical_water_flux.txt
(null)_0_water_balance.txt
(null)_0_SWC_begin.txt
(null)_0_SWC_end.txt
```

Depending on the code version, additional diagnostic outputs may include:

```text
debug_hydraulics_by_iter.txt
clamp diagnostics
SWC/SWP layer diagnostics
actual water uptake diagnostics
```

---

## 7. Key variables to inspect

### 7.1 Interface-level diagnostics

From the vertical water flux output:

```text
mean_flux_layersX_Y
mean_abs_flux_layersX_Y
mean_delta_swp_layersX_Y
mean_ks_interface_layersX_Y
net_volumetric_change_layersX_Y
gross_volumetric_change_layersX_Y
```

If the output still says `mean_ks_harmonic_layersX_Y`, check whether the underlying code actually uses harmonic or geometric mean.

### 7.2 Layer-level diagnostics

```text
mean_swp_layerX
mean_ks_layerX
SWC_layerX
theta_w_layerX
distance_to_Min_SWC
distance_to_Max_SWC
```

### 7.3 Water balance diagnostics

```text
precipitation
interception
throughfall
infiltration
runoff
leak
evaporation
transpiration by layer
actual water uptake
unmet transpiration
```

It is important to distinguish transpiration demand from actual water removed from the soil.

---

## 8. Interpretation criteria

### 8.1 SWP stability

Check whether SWP reaches extreme values such as:

```text
-1000 MPa
-2000 MPa
```

Such values are probably not useful ecologically and may indicate numerical artifacts, excessive gradients, or an overly conductive dry soil caused by high `Ks_floor`.

### 8.2 Vertical propagation

Check whether water exchange occurs only at interface `4_5`, or whether it propagates through:

```text
3_4
2_3
1_2
0_1
```

Possible interpretations:

- Strong `4_5` exchange but almost zero above: water table influences only the layer immediately above it.
- Strong exchange across multiple interfaces: water table or deep moisture influences upper soil layers.
- Very large upper-interface fluxes in `redprec`: may indicate excessive conductivity in dry layers.

### 8.3 Directionality

Use:

```text
directionality = abs(net_volumetric_change) / gross_volumetric_change
```

Interpretation:

- `directionality ~ 1`: mostly one-way flux.
- `directionality ~ 0`: strong back-and-forth exchange or oscillation.
- Large gross flux with low directionality may indicate numerical oscillation.

### 8.4 Climate sensitivity

A plausible result should show some distinction between `regclim` and `redprec`, especially in deepWT conditions.

Expected:

- lower SWC/SWP in upper layers under `redprec`;
- lower actual transpiration under `redprec` if water stress is active;
- stronger buffering in shallowWT than deepWT;
- different responses between sandy and clayey soils.

If `regclim` and `redprec` remain nearly identical, check:

1. whether precipitation reduction is actually applied;
2. whether throughfall and infiltration differ;
3. whether clamps or WT treatment are masking the difference;
4. whether actual water uptake differs from transpiration demand;
5. whether the simulation is still in spin-up.

---

## 9. Suggested R summaries

### 9.1 Flux summary by interface

```r
id_col <- intersect(
  c("model_name", "scenario", "experiment_name", "run_name", "scenario_base"),
  names(df_interface)
)[1]

flux_summary <- df_interface %>%
  filter(metric %in% c("gross_volumetric_change", "net_volumetric_change")) %>%
  group_by(.data[[id_col]], interface, metric) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    q95 = quantile(value, 0.95, na.rm = TRUE),
    max = max(value, na.rm = TRUE),
    .groups = "drop"
  )

flux_summary
```

### 9.2 Directionality

```r
directionality_summary <- df_interface %>%
  filter(metric %in% c("gross_volumetric_change", "net_volumetric_change")) %>%
  select(all_of(id_col), sim_year, interface, metric, value) %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  mutate(
    directionality = abs(net_volumetric_change) / gross_volumetric_change
  ) %>%
  group_by(.data[[id_col]], interface) %>%
  summarise(
    mean_gross = mean(gross_volumetric_change, na.rm = TRUE),
    median_gross = median(gross_volumetric_change, na.rm = TRUE),
    q95_gross = quantile(gross_volumetric_change, 0.95, na.rm = TRUE),
    max_gross = max(gross_volumetric_change, na.rm = TRUE),
    mean_net = mean(net_volumetric_change, na.rm = TRUE),
    median_directionality = median(directionality, na.rm = TRUE),
    .groups = "drop"
  )

directionality_summary
```

### 9.3 Log-scale plot for small fluxes

```r
df_interface %>%
  filter(metric == "gross_volumetric_change") %>%
  mutate(value_plot = pmax(value, 1e-30)) %>%
  ggplot(aes(sim_year, value_plot)) +
  geom_line(linewidth = 0.3) +
  facet_grid(.data[[id_col]] ~ interface, scales = "free_y") +
  scale_y_log10() +
  theme_bw()
```

This plot helps determine whether upper-interface fluxes are truly zero or only much smaller than the `4_5` flux.

---

## 10. Notes on current interpretation

Current preliminary results suggest:

1. `theta_w_min` is not the main driver of the observed behavior.
2. `Ks_floor` strongly controls vertical connectivity.
3. `Ks_floor = 1e-10` may make dry layers too conductive.
4. `Ks_floor = 1e-20` may make the profile too disconnected above the water table.
5. `Ks_floor = 1e-16` still appears strongly disconnected above interface `4_5` in deepWT sandy conditions.
6. Intermediate values such as `1e-12` and `1e-14` should be tested.
7. The deepWT sandy case should be treated as the primary stress test for hydraulic locking.
8. If the model remains extremely sensitive to interface conductivity, vertical layer discretization and timestep length should be revisited.

---

## 11. Open issues

The following issues still need to be checked before interpreting model outputs ecologically:

- Is precipitation reduction applied only through climate input files, or is there still a hard-coded precipitation multiplier in the C++ code?
- Is the water table treated as both a source and a sink, or only as a source?
- Are saturated layers below the WT allowed to receive recharge?
- Are full thick layers being saturated when only part of the layer lies below the WT?
- Are SWC clamps creating or destroying non-negligible amounts of water?
- Is transpiration output reporting demand or actual water removed from soil?
- Are differences between `regclim` and `redprec` visible in throughfall, infiltration, SWC, SWP, and actual transpiration?
- Is the simulation long enough to separate spin-up from equilibrium behavior?

---

## 12. Recommended next step

Run the following controlled tests first:

```text
deepWT_sandy
veg
regclim and redprec
theta_w_min = 1e-3
interface K = geometric
Ks_floor = 1e-12, 1e-14, 1e-16
```

Then compare:

```text
SWP by layer
gross and net volumetric change by interface
directionality
mean interface conductivity
actual water uptake
water balance closure
```

A good candidate `Ks_floor` should avoid both extremes:

- not so high that dry layers remain artificially conductive;
- not so low that the entire profile above `4_5` becomes completely disconnected;
- able to produce plausible differences between `regclim` and `redprec`;
- stable over long simulations without extreme SWP values or large artificial water creation/destruction.
