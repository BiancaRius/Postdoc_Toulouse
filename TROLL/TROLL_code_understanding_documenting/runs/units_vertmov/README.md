# Testing unit-consistency corrections in the capillary rise / vertical water flux routine

## Purpose

This set of simulations is intended to test recent corrections made to the capillary rise / vertical water flux routine in TROLL.

The main goal is to check whether previous water mass-balance problems were caused, at least partly, by inconsistent treatment of units inside the capillary rise function. In particular, the previous implementation mixed:

- volumetric soil water content (`SWC3D`, m³ water / m³ soil),
- water height or water depth (`water_height_upward`, m water),
- absolute water volume (`water_change_vol`, m³ water).

This could lead to incorrect updates of `SWC3D`, artificial overshoots beyond physical soil water limits, and subsequent creation or destruction of water by the numerical clamp.

## Main issue being tested

The capillary rise / vertical water flux routine computes water movement between adjacent soil layers. However, the state variable updated by the model is `SWC3D`, which is a volumetric water content.

Therefore, any water transfer computed as an absolute water volume must be converted back to a volumetric water content change before updating `SWC3D`.

The correct conversion is:

```cpp
delta_swc = water_change_vol[l] / layer_volume[l];
SWC3D[l][d] += delta_swc;