# Reserva Ducke-like soil and water table experiments

## Objective

This set of simulations aims to mimic the soil and water table configuration observed in Reserva Ducke, Manaus, where soil texture and water table depth are strongly linked to topography.

In this landscape, shallow water table conditions are associated with sandy soils, mainly due to erosion and redistribution processes from upper land soils. In contrast, deeper water table conditions are associated with clayey soils, typically found in higher topographic positions.

These runs are designed to test how this topographic coupling between soil texture and water table depth affects soil water dynamics, vegetation functioning, and functional trait distributions in TROLL.

## Experimental design

The initial simulations will be run under regular climate conditions.

The main configurations are:

- **Shallow water table + sandy soil**
  - Represents lower topographic positions.
  - Mimics areas where sandy soils are present due to erosion processes from upper land soils.
  - Expected to show stronger interaction between the rooting zone and groundwater.

- **Deep water table + clayey soil**
  - Represents upper topographic positions.
  - Mimics areas with deeper groundwater and more clay-rich soils.
  - Expected to show weaker direct groundwater influence on the upper soil layers.

## Initial climate forcing

The first block of simulations will use the regular climate scenario.

Reduced precipitation or other climate perturbation scenarios may be added later, after the baseline behavior of the coupled soil–water table configurations is understood.

## Analysis plan

The simulations will be analyzed from three complementary perspectives.

### 1. Soil water dynamics

The hydrological analysis will focus on how soil texture and water table depth control water redistribution through the soil profile.

Key variables to inspect include:

- Soil water content by layer
- Soil water potential by layer
- Vertical water fluxes between layers
- Gross and net volumetric water changes by interface
- Hydraulic conductivity by interface
- Water table influence on upper soil moisture
- Evidence of upward flux or hydraulic disconnection

### 2. Biogeochemical variables

The second block of analysis will evaluate whether differences in soil and water table configuration affect ecosystem-level functioning.

Variables of interest include:

- Aboveground biomass
- Basal area
- Gross primary productivity
- Net primary productivity
- Respiration components
- Litterfall
- Transpiration
- Water stress effects on vegetation dynamics

### 3. Trait distribution

The third block of analysis will evaluate whether hydrological conditions filter plant strategies differently across the simulated topographic settings.

Trait-related outputs may include:

- Trait distributions through time
- Shifts in dominant strategies
- Changes in functional composition
- Potential filtering by water availability, rooting depth, or drought tolerance
- Differences between shallow WT + sandy soil and deep WT + clayey soil communities

