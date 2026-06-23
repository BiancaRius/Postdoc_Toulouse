# Long-term comparison: original TROLL code vs implementation with deep water table

## Objective

The goal of this analysis is to compare long-term simulations between the original TROLL code and the modified version of the code including the new hydrological implementations.

The comparison will be performed using 1000-year simulations under regular climate conditions, using the Paracou climate forcing and sandy soil.

The main comparison is:

**Original code vs deep water table implementation**

## Simulation setup

The simulations are based on the following configuration:

* Simulation length: 1000 years
* Climate: regular climate
* Climate forcing: Paracou
* Soil type: sandy soil
* Comparison:

  * Original TROLL code
  * Modified TROLL code with deep water table implementation

## Output directories

### Original code

```bash
/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/original_code_TROLL/original_code_regclim_longterm/output
```

### Deep water table implementation

```bash
/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/TROLL_code_understanding_documenting/runs/longterm_orig_code_vs_wt/longterm_deep_sandy_regclim/output
```

## Purpose of the comparison

This comparison is intended to evaluate whether the new implementation produces plausible long-term ecosystem dynamics when compared to the original TROLL code.

In particular, the analysis will focus on checking whether the simulations:

* remain numerically stable over long timescales;
* produce realistic values for key ecosystem variables;
* reach a stable or quasi-stable state after the initial transient period;
* show consistent long-term behavior between the original and modified versions of the code.

## Variables to compare

The main variables to be compared may include:

* Aboveground biomass
* Basal area
* Gross primary productivity
* Net primary productivity
* Stem density
* Litterfall
* Evaporation
* Transpiration
* Soil water content
* Soil water potential
* Water balance components

## Expected output

The comparison will be based on time series plots and summary statistics of the long-term simulations.

The goal is not necessarily to obtain identical results between the original and modified versions, but to verify whether the implementation behaves consistently and produces ecologically plausible outputs under the same climate and soil conditions.

## Notes

This analysis is part of the validation process of the hydrological implementation in TROLL. The first step is to compare the original code with the deep water table scenario under regular climate and sandy soil conditions before exploring additional scenarios.
