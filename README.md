# Multiple Tidal Lagoons Hydrodynamics Model — UK

A 2D hydrodynamic model for simulating multiple idealised tidal range energy conversion (TREC) lagoons in the Severn Estuary and wider West UK coastal region, built on the [Thetis](https://thetisproject.org/) finite element ocean model.

---

## Scientific Background

The Severn Estuary has one of the largest tidal ranges in the world (up to ~14 m at spring tide), making it a prime candidate for tidal range energy extraction. Tidal lagoons — partially enclosed coastal impoundments fitted with low-head hydro turbines — can exploit this resource by generating electricity during both the ebb and flood phases of the tide.

While individual lagoon designs have been studied extensively, the **cumulative hydrodynamic impact** of deploying multiple lagoons simultaneously remains poorly understood. Modifying the tidal prism in one lagoon inevitably perturbs water levels and currents throughout the estuary, potentially affecting the performance of neighbouring lagoons and the broader coastal environment.

This model addresses that gap by simulating up to **seven idealised TREC lagoons** in the West UK simultaneously, allowing researchers to:

- Quantify inter-lagoon hydrodynamic interactions.
- Assess changes to tidal range, currents, and residual circulation.
- Evaluate aggregate power generation and capacity factors.
- Inform policy decisions on consenting multiple tidal range projects.

The work follows the consistent design approach described in Mackie *et al.* — *"Assessing impacts of tidal power lagoons of a consistent design"* — where each lagoon is parameterised identically to isolate location-dependent effects.
