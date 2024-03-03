# # Springpot
#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/model_elements.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/model_elements.ipynb)

using RHEOS
## include a helper function for plotting
include("assets/plothelper.jl");
#-

# By typing the name of the model, it is possible to visualise its graphical representation and its parameters.

Springpot

# #### Constitutive Equation
# ```math
# \sigma(t) = c_{\beta} \frac{d^\beta \epsilon(t)}{dt^\beta}
# ```

# ```math
# \text{for}\; \ 0 \leq \beta \leq 1
# ```

# #### Relaxation Modulus
# ```math
# G(t) = \frac{c_{\beta} }{\Gamma(1-\beta)} t^{-\beta}
# ```

# #### Creep Modulus
# ```math
# J(t) = \frac{1}{c_\beta \Gamma(1+\beta)}t^\beta
# ```

# #### Storage Modulus
# ```math
# G^{\prime}(\omega) = c_\beta \omega^\beta \cos(\frac{\pi}{2}\beta)
# ```

# #### Loss Modulus
# ```math
# G^{\prime\prime}(\omega) = c_\beta \omega^\beta \sin(\frac{\pi}{2}\beta)
# ```

# # Spring

# When β = 0 the springpot specialises to a spring.

Spring

# # Dashpot

# When β = 1 the springpot specialises to a dashpot.

Dashpot

# ## Qualitative Behaviours of the Moduli

models = Vector{RheoModel}()

## Spring
push!(models, RheoModel(Spring, k = 1.0))

## plot moduli for varying β
for beta in [0.2, 0.5, 0.8]
    
    push!(models, RheoModel(Springpot, cᵦ = 1.0, β = beta))
 
end

## Dashpot
push!(models, RheoModel(Dashpot, η = 1.0))

plotmodel(models, ymaxG = 2.0)
