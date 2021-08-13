# Architecture

```@meta
EditURL = "https://github.com/JuliaRheology/RHEOS.jl/edit/master/docs/src/architecture.md"
```

RHEOS is built around four main data types:
+ [`RheoTimeData`](@ref): Contains time (t), stress (σ), strain (ϵ) data
+ [`RheoFreqData`](@ref): Contains frequency (ω), storage (Gp) and loss (Gpp) moduli
+ [`RheoModelClass`](@ref): Contains the model's name, parameters, and the expressions of relaxation (G), creep (J), storage (Gp) and loss (Gpp) moduli as functions of the model parameters
+ [`RheoModel`](@ref): Similar to `RheoModelClass`, but actual numbers are substituted (hard-coded) into the expressions for relaxation, creep, storage and loss moduli.

A common RHEOS workflow is illustrated in the figure below. Experimental time-domain viscoelastic data is fitted to a viscoelastic model. This model (with the fitted parameters) is then used to make a prediction of the behaviour under the fitted loading conditions so that its accuracy can be qualitatively assessed. Similarly, the fitted model can be used to simulate the behaviour of the same material under any loading conditions (different from the fitted ones).

![High level schematic of a fitting and prediction workflow from experimental data.](./assets/diagram.svg)

So, having seen the figure above we can now think about our rheology data analysis workflow with reference to RHEOS' types and functions. Experimental data are imported into a `RheoTimeData` struct (or `RheoFreqData` for dynamic experimental data). This is then fitted to a `RheoModelClass`, which is the model struct used when the model's parameters are not fixed. The output of this process is a `RheoModel` (which is a `RheoModelClass` in which the parameters have been fixed to specific values). In the prediction step, the fitted `RheoModel` is combined with partial data (`RheoTimeData` with time and either only stress or strain) and fills in the missing data column, resulting in a complete data set (complete `RheoTimeData`).
