# Sensitivity tests for Ks_floor

This folder contains a set of TROLL simulations designed to test the effect of different minimum hydraulic conductivity values (`Ks_floor`) on vertical soil water redistribution.

## Goal

The goal of these tests is to evaluate whether changing `Ks_floor` affects:

- hydraulic locking;
- upward water movement from the water table;
- soil water content dynamics;
- soil water potential dynamics;
- numerical stability of the model;
- water balance consistency.

## Tested values

The tested `Ks_floor` values are:

```text
1e-10
1e-12
1e-14
1e-16
1e-18
1e-20