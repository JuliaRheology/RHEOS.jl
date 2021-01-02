# # Fractional Maxwell
#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/model_maxwell.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/model_maxwell.ipynb)

using RHEOS
## include a helper function for plotting
include("assets/plothelper.jl");
#-

Fract_Maxwell

# #### Constitutive Equation
# ```math
# \sigma(t) + \frac{c_\alpha}{c_\beta} \frac{d^{\alpha-\beta} \sigma(t)}{dt^{\alpha-\beta}}= c_{\alpha} \frac{d^\alpha \epsilon(t)}{dt^\alpha}
# ```

# ```math
# \text{for}\; \ 0 \leq \beta \leq \alpha \leq 1
# ```

# #### Relaxation Modulus
# ```math
# G(t) = c_\beta t^{-\beta} E_{\alpha-\beta,1-\beta}\left(-\frac{c_\beta}{c_\alpha} t^{\alpha-\beta}\right)
# ```

# #### Creep Modulus
# ```math
# J(t) = \frac{1}{c_\alpha \Gamma(1+\alpha)}t^\alpha+\frac{1}{c_\beta \Gamma(1+\beta)}t^\beta
# ```

# #### Storage Modulus
# ```math
# G^{\prime}(\omega) = \frac{\left(c_\beta \omega^\beta\right)^2 \cdot c_\alpha \omega^\alpha \cos(\alpha \frac{\pi}{2}) + \left(c_\alpha \omega^\alpha\right)^2 \cdot c_\beta \omega^\beta \cos(\beta \frac{\pi}{2})}{\left(c_\alpha \omega^\alpha\right)^2+\left(c_\beta \omega^\beta\right)^2+2c_\alpha \omega^\alpha \cdot c_\beta \omega^\beta \cos((\alpha-\beta)\frac{\pi}{2})}
# ```

# #### Loss Modulus
# ```math
# G^{\prime\prime}(\omega) = \frac{\left(c_\beta \omega^\beta\right)^2 \cdot c_\alpha \omega^\alpha \sin(\alpha \frac{\pi}{2}) + \left(c_\alpha \omega^\alpha\right)^2 \cdot c_\beta \omega^\beta \sin(\beta \frac{\pi}{2})}{\left(c_\alpha \omega^\alpha\right)^2+\left(c_\beta \omega^\beta\right)^2+2c_\alpha \omega^\alpha \cdot c_\beta \omega^\beta \cos((\alpha-\beta)\frac{\pi}{2})}
# ```

# # Fractional (Spring) Maxwell

FractS_Maxwell
#-

models = Vector{RheoModel}()

## plot moduli for varying α
for alpha in [0.1, 0.3, 0.5, 0.7, 0.9]
    
    push!(models, RheoModel(FractS_Maxwell, (cₐ = 1.0, a = alpha, k = 1.0)))
 
end

plotmodel(models)

# # Fraction (Dashpot) Maxwell

FractD_Maxwell
#-

models = Vector{RheoModel}()

## plot moduli for varying β
for beta in [0.1, 0.3, 0.5, 0.7, 0.9]
    
    push!(models, RheoModel(FractD_Maxwell, (η = 10, cᵦ= 1.0, β = beta)))
 
end

plotmodel(models, ymaxG = 2.0)

# # Maxwell

Maxwell
#-

models = Vector{RheoModel}()

## plot moduli for varying k
for k in [5.0, 10.0, 20.0]
    
    push!(models, RheoModel(Maxwell, (η = 10, k = k)))
 
end

plotmodel(models)