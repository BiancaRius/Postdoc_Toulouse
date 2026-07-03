# Debug runs after unit corrections in vertical water movement

## Why these runs are being performed

These runs are being performed because several unit inconsistencies were identified in the implementation of vertical water movement and capillary rise in TROLL.

The objective of these simulations is not yet to interpret long-term ecological or hydrological results, but to verify whether the corrected implementation is now numerically stable, dimensionally consistent, and physically interpretable.

Previous runs should therefore be treated as debugging runs rather than final model results.

## Unit inconsistencies identified

### 1. Hydraulic conductivity units

The saturated hydraulic conductivity `Ksat` is calculated from the Cosby/Marthews formulation in:

```text
mm/s
```

However, the Darcy flux calculation requires hydraulic conductivity in:

```text
m/s
```

Previously, `Ks_cap` and the interface conductivity were being used in Darcy’s law as if they were already in `m/s`.

Corrected convention:

```text
Ksat originally calculated = mm/s
Ks_cap after conversion    = m/s
interface Ks               = m/s
q_cap                      = m/s
```

Correction applied:

```cpp
Ks_cap[l][d] = Ks_cap[l][d] * 1e-3f; // mm/s -> m/s
```

The water table case also needs this conversion:

```cpp
Ks_cap[l][d] = Ksat[l] * 1e-3f; // mm/s -> m/s
```

### 2. Ks floor interpretation

Because `Ks_cap` is now converted to `m/s` before applying the lower bound, the `ksfloor` value is interpreted in:

```text
m/s
```

Therefore:

```cpp
if (Ks_cap[l][d] < 1e-14f) {
    Ks_cap[l][d] = 1e-14f; // m/s
}
```

If the floor were applied before conversion, the equivalent value would be different:

```text
1e-14 m/s = 1e-11 mm/s
```

The corrected implementation applies the floor after conversion to `m/s`.

### 3. Soil water storage units

The following variables are stored as water volumes:

```text
SWC3D
Min_SWC
Max_SWC
FC_SWC
water_change_vol
water_change_cap
```

Unit:

```text
m3 water
```

Previously, some parts of the code treated these variables as if they were volumetric water contents in `m3/m3`.

Corrected convention:

```text
SWC3D = m3 water
Min_SWC = m3 water
Max_SWC = m3 water
FC_SWC = m3 water
```

Therefore, water volume changes should be applied directly:

```cpp
SWC3D[l][d] += water_change_vol[l]; // m3 water
```

and not as:

```cpp
SWC3D[l][d] += water_change_vol[l] / layer_volume[l];
```

### 4. Receiver and donor capacities

The variables:

```text
max_gain
max_loss
receiv_capacity
donor_capacity
```

are also water volumes.

Unit:

```text
m3 water
```

Because `Max_SWC`, `FC_SWC`, `Min_SWC`, and `SWC3D` are already in `m3 water`, the differences:

```cpp
Max_SWC[l] - SWC3D[l][d]
SWC3D[l][d] - Min_SWC[l]
```

are already in:

```text
m3 water
```

Previously, these values were multiplied again by `layer_volume`, which produced an inconsistent unit:

```text
m3 water * m3 soil = m6
```

Corrected convention:

```cpp
receiv_capacity[l] = max_gain[l]; // m3 water
donor_capacity[l]  = max_loss[l]; // m3 water
```

### 5. Infiltration into the top layer

In the unified vertical flux scheme, throughfall input `in` is available in:

```text
m3 water
```

The amount of water layer 0 can receive is:

```cpp
receiv_l0 = Max_SWC[0] - SWC3D[0][d];
```

Unit:

```text
m3 water
```

Previously, this receiving capacity was multiplied by `layer_volume`, and the infiltrated volume was divided by `layer_volume` before updating `SWC3D`.

Corrected convention:

```cpp
float receiv_l0 = std::max(0.0f, Max_SWC[0] - SWC3D[0][d]); // m3 water
float infil0_vol = fminf(in, receiv_l0);                    // m3 water

SWC3D[0][d] += infil0_vol;                                  // m3 water
```

No multiplication or division by `layer_volume` is needed.

### 6. Potential transferred water volume

The Darcy flux produces:

```text
q_cap = m/s
```

The water height moved during one timestep is:

```cpp
water_height_upward[l][d] = q_cap[l][d] * delta_t_sec;
```

Unit:

```text
m
```

The potential transferred volume is then:

```cpp
vol_potential_transfer = water_height_upward[l][d] * voxel_area;
```

Unit:

```text
m * m2 = m3 water
```

This volume can then be limited by donor and receiver capacities, which are also in `m3 water`.

### 7. Clamp diagnostic conversion

The clamp diagnostic computes:

```cpp
difference = SWC3D_after_clamp - SWC3D_before_clamp;
```

Since `SWC3D` is in `m3 water`, then:

```text
difference = m3 water
```

Previously, the equivalent water depth was calculated as if `difference` were a volumetric water content change in `m3/m3`.

Corrected conversion:

```cpp
delta_water_mm = difference / voxel_area * 1000.0;
```

because:

```text
m3 water / m2 = m water
m water * 1000 = mm water
```

### 8. Clamp and mass balance

The safety clamp improves numerical stability by preventing `SWC3D` from going outside physical limits. However, any clamp that changes `SWC3D` can create or remove a small amount of water.

Therefore, clamp-created and clamp-destroyed water are now diagnosed explicitly.

Relevant diagnostics:

```text
annual_clamp_created_swc
annual_clamp_destroyed_swc
annual_clamp_created_mm
annual_clamp_destroyed_mm
```

These values should remain negligible. If they are not negligible, the clamp is affecting the water balance and must be reviewed.

## First debug run

The first corrected run should use a simple scenario where vertical water movement is expected to be easy to detect:

```text
shallow water table + sandy soil + regular climate
```

Recommended scenario name:

```text
shallow_sandy_regclim
```

Recommended duration:

```text
5 years
```

If the model runs without numerical problems, the same scenario can be extended to 10 years.

## Purpose of the first debug run

This first run is intended to check:

```text
1. whether the model runs without crashing
2. whether SWC3D remains within physical limits
3. whether layers below the water table remain saturated
4. whether q_cap has plausible values in m/s
5. whether water_height_upward has plausible values in m
6. whether water_upward_vol and water_change_vol are in m3
7. whether donor and receiver capacities are respected
8. whether clamp-created or clamp-destroyed water is negligible
9. whether artificial hydraulic locking is reduced after the unit corrections
```

## Suggested folder name

```text
debug_vertical_flux_units_shallow_sandy_regclim_5yr
```
