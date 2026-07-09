Water Table Depth Tests with Fixed Soil Profile
Objective

This set of simulations tests the effect of water table depth (WTD) on soil water dynamics, vertical water movement, capillary rise, and vegetation responses in TROLL.

The main goal is to avoid confounding water table depth with changes in soil profile structure. Therefore, all simulations use the same six-layer soil profile. The only factor changing among these tests is the imposed water table depth.

Soil profile

The model uses six soil layers with the following thicknesses:

layer_thickness = 0.1, 0.5, 1.0, 2.4, 6.0, 10.0 m

The code computes the cumulative depth of each layer as:

cumulative_depth += layer_thickness[l];
layer_depth[l] = cumulative_depth;

Therefore, layer_depth[l] corresponds to the lower boundary of each soil layer.

The resulting soil profile is:

Layer	Thickness (m)	Top depth (m)	Bottom depth / layer_depth[l] (m)
0	0.1	0.0	0.1
1	0.5	0.1	0.6
2	1.0	0.6	1.6
3	2.4	1.6	4.0
4	6.0	4.0	10.0
5	10.0	10.0	20.0

The total soil profile depth is therefore 20 m.

Water table depth definition

Water table depth is interpreted as the upper boundary of the saturated zone.

In the current implementation, a layer is treated as saturated when its lower boundary is deeper than the imposed water table depth:

if (_WATER_TABLE == 1 && layer_depth[l] > WTD) {
    soil_phi3D_cap[l][d] = 0.0f;
    soil_phi3D[l][d] = soil_phi3D_cap[l][d];
    theta_w_cap = 1.0f;
    Ks_cap[l][d] = Ksat[l] * 1e-3f;
    continue;
}

Because layer_depth[l] is the cumulative lower boundary of the layer, the selected WTD values are chosen to correspond exactly to the top of the first saturated layer.

Water table scenarios

The planned WTD scenarios are:

WTD = 10.0, 4.0, 1.6, 0.6, 0.0 m

These values correspond to saturation starting at progressively shallower layers.

Scenario label	WTD (m)	Saturated layers	Interpretation
WT10 / deepWT	10.0	layer 5	Deep water table baseline
WT4	4.0	layers 4 and 5	Saturated zone starts at 4.0 m
WT1p6	1.6	layers 3, 4, and 5	Saturated zone starts at 1.6 m
WT0p6	0.6	layers 2, 3, 4, and 5	Saturated zone starts at 0.6 m
WT0	0.0	layers 0, 1, 2, 3, 4, and 5	Fully saturated soil profile

The baseline scenario is:

WTD = 10.0 m

In this baseline, only the deepest layer, from 10 to 20 m, is saturated.

Rationale

The previous comparison between shallow and deep water table scenarios was potentially confounded by differences in effective soil layer thickness and saturated-layer configuration.

In this new setup:

The soil profile is fixed.
The layer thicknesses are fixed.
The total soil depth is fixed.
Only WTD changes among scenarios.

This allows differences among simulations to be attributed more directly to water table depth.