# Water Table Depth Tests with Fixed Soil Profile

## Objective

This set of simulations tests the effect of water table depth (WTD) on soil water dynamics, vertical water movement, capillary rise, and vegetation responses in TROLL.

The main goal is to avoid confounding water table depth with changes in soil profile structure. Therefore, all simulations use the same six-layer soil profile. The only factor changing among scenarios is the imposed water table depth.

---

## Soil profile

The model uses six soil layers with the following thicknesses:

```text
layer_thickness = 0.1, 0.5, 1.0, 2.4, 6.0, 10.0 m
```

In the code, cumulative layer depth is calculated as:

```cpp
cumulative_depth += layer_thickness[l];
layer_depth[l] = cumulative_depth;
```

Therefore, `layer_depth[l]` corresponds to the lower boundary of each soil layer.

The resulting soil profile is:

| Layer | Thickness (m) | Top depth (m) | Bottom depth / `layer_depth[l]` (m) |
| ----- | ------------: | ------------: | ----------------------------------: |
| 0     |           0.1 |           0.0 |                                 0.1 |
| 1     |           0.5 |           0.1 |                                 0.6 |
| 2     |           1.0 |           0.6 |                                 1.6 |
| 3     |           2.4 |           1.6 |                                 4.0 |
| 4     |           6.0 |           4.0 |                                10.0 |
| 5     |          10.0 |          10.0 |                                20.0 |

The total soil profile depth is therefore 20 m.

---

## Water table depth definition

Water table depth is interpreted as the upper boundary of the saturated zone.

In the current implementation, a layer is treated as saturated when its lower boundary is deeper than the imposed WTD:

```cpp
if (_WATER_TABLE == 1 && layer_depth[l] > WTD) {
    soil_phi3D_cap[l][d] = 0.0f;
    soil_phi3D[l][d] = soil_phi3D_cap[l][d];
    theta_w_cap = 1.0f;
    Ks_cap[l][d] = Ksat[l] * 1e-3f;
    continue;
}
```

Because the condition uses:

```cpp
layer_depth[l] > WTD
```

the selected WTD values are placed at the top of the first saturated layer.

---

## Water table scenarios

The planned WTD scenarios are:

```text
WTD = 10.0, 4.0, 1.6, 0.6, 0.1, 0.0 m
```

These values correspond to saturation starting at progressively shallower layers.

| Scenario label    | WTD (m) | Saturated layers            | Interpretation                  |
| ----------------- | ------: | --------------------------- | ------------------------------- |
| `WT10` / `deepWT` |    10.0 | layer 5                     | Saturated zone starts at 10.0 m |
| `WT4`             |     4.0 | layers 4 and 5              | Saturated zone starts at 4.0 m  |
| `WT1p6`           |     1.6 | layers 3, 4, and 5          | Saturated zone starts at 1.6 m  |
| `WT0p6`           |     0.6 | layers 2, 3, 4, and 5       | Saturated zone starts at 0.6 m  |
| `WT0p1`           |     0.1 | layers 1, 2, 3, 4, and 5    | Saturated zone starts at 0.1 m  |
| `WT0`             |     0.0 | layers 0, 1, 2, 3, 4, and 5 | Fully saturated soil profile    |

The baseline scenario is:

```text
WTD = 10.0 m
```

In this baseline, only the deepest layer, from 10 to 20 m, is saturated.

---

## Expected saturated layers

Given the current rule:

```cpp
layer_depth[l] > WTD
```

the expected saturated layers are:

```text
WTD = 10.0 -> layer 5
WTD = 4.0  -> layers 4 and 5
WTD = 1.6  -> layers 3, 4, and 5
WTD = 0.6  -> layers 2, 3, 4, and 5
WTD = 0.1  -> layers 1, 2, 3, 4, and 5
WTD = 0.0  -> layers 0, 1, 2, 3, 4, and 5
```

This should be verified in the model output or diagnostic prints before running the full experiment.

---

## Rationale

Previous comparisons between shallow and deep water table scenarios were potentially confounded by differences in effective soil layer thickness and saturated-layer configuration.

In this new setup:

```text
The soil profile is fixed.
The layer thicknesses are fixed.
The total soil depth is fixed.
Only WTD changes among scenarios.
```

This allows differences among simulations to be attributed more directly to water table depth.

The main question is:

```text
How does progressively shallower water table depth affect soil water availability,
vertical hydraulic connectivity, capillary rise, transpiration, and forest dynamics?
```

---

## Important implementation note

The current implementation assumes that WTD coincides with a layer boundary.

The selected WTD values:

```text
10.0, 4.0, 1.6, 0.6, 0.1, 0.0 m
```

correspond exactly to the top of layers 5, 4, 3, 2, 1, and 0, respectively.

Intermediate WTD values should be avoided unless partial saturation within a layer is explicitly implemented.

For example, with the current logic, if an intermediate WTD such as 2.0 m were used, the entire layer 3 would be treated as saturated because its lower boundary is deeper than 2.0 m. This would shift the effective saturated zone upward to the top of layer 3.


```


---

## Main diagnostics to analyze

### 1. Soil water content

Variables to check:

```text
SWC by layer
theta_w by layer
distance to Min_SWC
distance to Max_SWC
```

Expected pattern:

* saturated layers should remain at or near saturation;
* shallower WTD should increase water availability in upper and intermediate layers;
* `WT10` should represent the deepest baseline, with only layer 5 saturated.

---

### 2. Soil water potential

Variables to check:

```text
SWP by layer
soil_phi3D by layer
soil_phi3D_cap by layer
hydraulic head by layer
delta hydraulic head by interface
```

Expected pattern:

* saturated layers should have matric potential set to zero;
* unsaturated layers should respond to rainfall, drainage, root uptake, and capillary rise;
* shallower WTD should reduce water stress in layers closer to the root zone.

---

### 3. Hydraulic conductivity

Variables to check:

```text
Ks_cap by layer
Ks_harmonic by interface
```

Expected pattern:

* saturated layers should have `Ks_cap = Ksat * 1e-3`;
* unsaturated layers should have conductivity controlled by soil moisture;
* deep WTD may show low connectivity if the unsaturated layer above the saturated zone becomes dry.

---

### 4. Vertical fluxes

Variables to check:

```text
mean_flux
mean_abs_flux
gross volumetric change
net volumetric change
directionality = abs(net) / gross
```

Questions to evaluate:

```text
Does upward flux increase as WTD becomes shallower?
Does WT10 show hydraulic disconnection after spin-up?
Is vertical flux limited by hydraulic gradient, conductivity, or donor/receiver capacity?
```

---

### 5. Transpiration by layer

Variables to check:

```text
Transpiration by layer
Total transpiration
Root distribution by layer
```

Expected pattern:

* shallower WTD may increase access to water;
* very shallow WTD may reduce root access if hypoxia/anoxia restrictions are implemented;
* deeper WTD may force trees to rely more strongly on water stored in upper layers after rainfall.

---

### 6. Ecosystem outputs

Variables to check:

```text
GPP
NPP
AGB
BA
mortality
recruitment
litterfall
```

Questions to evaluate:

```text
Does shallower WTD increase productivity by improving water access?
Does very shallow WTD reduce productivity because of saturated/anoxic soil?
Is there an intermediate WTD that maximizes productivity?
```

---

## Suggested first analysis

For the first round, compare only the regular climate scenarios after spin-up.

Recommended filter:

```text
year >= 101
```

Recommended plots:

1. Annual mean SWC by layer;
2. Annual mean SWP by layer;
3. Annual transpiration by layer;
4. Total annual transpiration;
5. Vertical flux by interface;
6. `Ks_harmonic` by interface;
7. Gross and net vertical water movement;
8. Directionality index: `abs(net) / gross`;
9. GPP, NPP, AGB, and BA through time;
10. Post-spinup boxplots or means across WTD scenarios.

---

## Expected qualitative behavior

The working expectations are:

1. `WT10` is the deep water table baseline, with only layer 5 saturated.
2. Moving the WTD upward should increase the potential for upward water supply.
3. Intermediate WTD may increase productivity by improving plant water access.
4. Very shallow WTD may reduce productivity if saturated soil restricts roots or causes hypoxia/anoxia.
5. Deep WTD may show hydraulic disconnection if the unsaturated layers above the saturated zone become too dry and hydraulic conductivity collapses.
6. Because the soil profile is fixed, differences among scenarios should mainly reflect water table depth rather than differences in soil discretization.

---

## Summary

This experiment uses a fixed six-layer soil profile:

```text
layer_thickness = 0.1, 0.5, 1.0, 2.4, 6.0, 10.0 m
```

which gives cumulative layer depths:

```text
layer_depth = 0.1, 0.6, 1.6, 4.0, 10.0, 20.0 m
```

The WTD scenarios are:

```text
10.0, 4.0, 1.6, 0.6, 0.1, and 0.0 m
```

Using the current rule:

```cpp
layer_depth[l] > WTD
```

these correspond to saturation starting respectively at layers:

```text
5, 4, 3, 2, 1, and 0
```

This setup allows a clean comparison of water table depth effects while keeping the soil profile constant across all simulations.
