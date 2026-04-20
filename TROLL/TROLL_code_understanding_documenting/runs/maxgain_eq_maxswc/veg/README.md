# Hydraulic tests for vertical water flux and capillary rise

## Overview

This set of runs was designed to test the recent implementation of vertical soil water fluxes in TROLL, with special attention to numerical stability and to warnings that have appeared in the log during simulations.

The main goal is to understand whether the current formulation produces physically coherent soil water dynamics under different combinations of:
- water table depth
- soil texture
- vegetation presence/absence
- precipitation reduction
- numerical thresholds used for hydraulic calculations

A particular focus of these tests is the treatment of warnings associated with the `CapillaryRise()` routine, especially messages indicating that soil water content became **out of bounds** after the capillary-rise update.

## Why these tests are needed

Recent simulations have produced repeated warnings such as:

> `Warning: layer 0 is above WTD but SWC3D out of bounds after capillary rise`

These warnings suggest that, after the vertical flux update, the soil water content of a layer above the water table may exceed the expected range used by the code.

At the moment, the implementation of vertical fluxes combines:
- explicit upward/downward transfer between layers
- special handling of layers below the water table
- bounds and safety checks applied after the update step

In the current code, the receiver capacity in `CapillaryRise()` is still computed from `FC_SWC` even when `_UNIFIED_VERT_WATER_FLUX == 1`, while layers below the water table are forced back to `Max_SWC` during the update, and warnings are printed when updated `SWC3D` falls outside the expected range. :contentReference[oaicite:0]{index=0} :contentReference[oaicite:1]{index=1}

Because of this, part of the present tests is not only to check for crashes or NaNs, but also to determine whether these warnings indicate:
1. a real physical inconsistency,
2. a numerical artifact,
3. or a mismatch between the physical limit that should be used in the unified Darcy scheme and the current diagnostic threshold.

## Main objectives

The runs in this folder were created to:

1. evaluate whether the new vertical water flux implementation behaves consistently across environmental settings;
2. identify under which conditions warnings become frequent;
3. distinguish harmless warnings from warnings associated with unstable model behavior;
4. test alternative numerical safeguards for very low soil moisture and conductivity values;
5. assess whether the current upper bound used for layers above the water table is appropriate in the unified vertical flux scheme.

## Key warnings under investigation

The main warnings of interest are those related to soil moisture bounds after capillary rise, especially:
- `SWC3D out of physical bounds immediately after update`
- `layer X is above WTD but SWC3D out of bounds after capillary rise`

These warnings are especially important because they may reveal that:
- a layer is receiving more water than allowed during one timestep,
- the receiver limit is too restrictive or inconsistent with the chosen scheme,
- or water is accumulating in a way that is numerically possible in the current code but not conceptually intended.

Special attention is given to the **out-of-boundary** warnings affecting layers above the water table, as these are the warnings that have appeared most clearly in recent tests.

## What is being compared across runs

The comparison among runs focuses on:
- shallow versus deep water table
- sandy versus clayey soil
- runs with and without vegetation
- standard versus reduced precipitation
- different threshold values for hydraulic variables such as `theta_w` and conductivity floors

## Variables and outputs to inspect

To interpret the warnings, the following outputs are especially relevant:
- `SWC3D`
- `soil water potential`
- layer-to-layer vertical fluxes
- net and gross changes in water volume between layers
- harmonic hydraulic conductivity at interfaces
- diagnostic counters for small `theta_w`, small `Ks`, and small harmonic `Ks`

When warnings appear, the most useful quantities to inspect are:
- the affected layer
- the current `SWC3D`
- its expected range
- the applied water change in the timestep
- receiver and donor capacities
- whether the warning occurs repeatedly in the same layer or propagates through time

## Working hypothesis

A current working hypothesis is that some of the warnings above the water table may be linked to the fact that the unified vertical flux formulation allows water redistribution more dynamically than the old bucket logic, while the diagnostic bound still relies on a stricter expectation for what an unsaturated layer should contain.

In that case, not every warning would necessarily imply a severe model failure. However, if warnings are followed by NaNs, unrealistic soil moisture profiles, or implausible plant water stress, they should be treated as a sign of instability and not only as a diagnostic message.

## Expected outcome of these tests

These runs are intended to help decide:
- whether the current implementation is robust enough to be kept as is,
- whether the warning logic should be revised,
- whether the upper limit for layers above the water table should be reformulated in the unified scheme,
- and whether additional numerical constraints are needed to avoid pathological soil water states.

## General note

The purpose of these tests is not only to suppress warnings, but to understand their origin and determine whether they reflect:
- real physical inconsistencies,
- purely numerical problems,
- or diagnostics that are currently too strict for the implemented formulation.