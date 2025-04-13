# Membrane Modeling Theory

## Overview
This document outlines the theory and methodology for modeling membrane separation processes, particularly for multicomponent systems. The approach involves discretizing the membrane, solving transport equations iteratively, and updating component permeabilities based on feed composition and temperature.

## Discretization and Permeability Calculation
The membrane is discretized into multiple steps using the following equation:

$$
dx = \frac{1}{n_{\text{points}}}, \quad dL = \frac{\text{length}}{n_{\text{points}}}
$$

Since permeability varies with feed composition but remains nearly constant at a fixed temperature (as supported by literature), the mole fraction of each component is used to calculate its permeability in Barrer. A Python function updates permeability values at each iteration for accuracy.

## Fugacity and Permeation Rate
The partial pressure and fugacity of each component are calculated using:

$$
f_{\text{feed}, i} = Z_{\text{feed}, i} \cdot P_{\text{feed}, i} \cdot 10^5
$$

$$
f_{\text{perm}, i} = Z_{\text{perm}, i} \cdot P_{\text{perm}, i} \cdot 10^5
$$

The driving force for permeation is the fugacity difference across the membrane. The permeation flow rate for each component is computed as:

$$
F_{\text{perm}, i} = \frac{P_i}{l} \cdot \max(f_{\text{feed}, i} - f_{\text{perm}, i}, 0) \cdot A \cdot dx
$$

## Total Flow Rates
The total permeate and retentate flow rates are obtained by summing individual component flow rates:

$$
F_{\text{ret}} = F_{\text{ret}} - \sum_{i} F_{\text{perm}, i}
$$

$$
F_{\text{perm}} = F_{\text{perm}} + \sum_{i} F_{\text{perm}, i}
$$

## Integration with Aspen Simulations
To avoid complexity, the Aspen model is not directly integrated into the Python simulation. Instead:
1. Feed data (composition, temperature, pressure) from Aspen is used as input for the Python membrane simulation.
2. Python simulation results (split fractions for the permeate stream) are fed into Aspen's separator (SEP) module.
3. This approach approximates the stream composition from Python within Aspen.

## Simulation Robustness and Optimization
The Python simulation was tested with varying feed compositions, revealing consistent split fractions. Key findings:
- **Optimal membrane area**: 103 m²  
- **Stage cut**: 0.9  
- **Hydrogen loss**: ~1.85%  

A contour graph illustrates the relationship between membrane area, stage cut, and hydrogen loss.

## References
- Pham et al. (2023) *Selective Permeability in Multicomponent Systems* (for permeability assumptions).
