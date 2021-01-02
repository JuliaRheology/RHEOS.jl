# # Fractional Poynting-Thomson
#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/notebooks/model_pt.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/model_pt.ipynb)

using RHEOS
## include a helper function for plotting
include("assets/plothelper.jl");
#-

Fract_PT
#-

# #### Constitutive Equation
# ```math
# \sigma(t) + \frac{c_\alpha}{c_\gamma} \frac{d^{\alpha-\gamma} \sigma(t)}{dt^{\alpha-\gamma}}+ \frac{c_\beta}{c_\gamma} \frac{d^{\beta-\gamma} \sigma(t)}{dt^{\beta-\gamma}}= c_{\alpha} \frac{d^\alpha \epsilon(t)}{dt^\alpha} + c_\beta \frac{d^\beta \epsilon(t)}{dt^\beta}
# ```

# ```math
# \text{for}\; \ 0 \leq \beta \leq \alpha \leq 1
# ```

# #### Relaxation Modulus
# ```math
# \tilde{G}(s) = \frac{1}{s}\frac{c_\gamma s^\gamma \cdot \left[c_\alpha s^\alpha + c_\beta s^{\beta}\right]}{c_\gamma s^\gamma+c_\alpha s^{\alpha}+c_\beta s^{\beta}}
# ```

# #### Creep Modulus
# ```math
# J(t)= \frac{t^{\alpha}}{c_\alpha} E_{\alpha-\beta,1+\alpha}\left(-\frac{c_\beta}{c_\alpha} t^{\alpha-\beta}\right) + \frac{1}{c_\gamma \Gamma(1+\gamma)}t^\gamma
# ```

# #### Storage Modulus
# ```math
# G^{\prime}(\omega) = \frac{c_\gamma \omega^\gamma \cos\left(\gamma \frac{\pi}{2}\right) \left[\left(c_\alpha \omega^\alpha\right)^2+\left(c_\beta \omega^\beta\right)^2 \right]+\left(c_\gamma \omega^\gamma\right)^2 \left[c_\alpha \omega^\alpha \cos\left(\alpha \frac{\pi}{2}\right)+c_\beta \omega^\beta \cos\left(\beta \frac{\pi}{2}\right) \right] + c_\alpha \omega^\alpha \cdot c_\beta\omega^\beta \cdot c_\gamma \omega^\gamma \left[\cos\left((\alpha-\beta-\gamma) \frac{\pi}{2}\right)+\cos\left((\beta-\alpha-\gamma) \frac{\pi}{2}\right) \right]}{\left(c_\alpha \omega^\alpha\right)^2+\left(c_\beta \omega^\beta\right)^2+\left(c_\gamma \omega^\gamma\right)^2+2c_\alpha \omega^\alpha \cdot c_\beta \omega^\beta \cos((\alpha-\beta)\frac{\pi}{2})+2c_\alpha \omega^\alpha \cdot c_\gamma \omega^\gamma \cos((\alpha-\gamma)\frac{\pi}{2})+2c_\beta \omega^\beta \cdot c_\gamma \omega^\gamma \cos((\beta-\gamma)\frac{\pi}{2})}
# ```

# #### Loss Modulus
# ```math
# G^{\prime\prime}(\omega) = \frac{c_\gamma \omega^\gamma \sin\left(\gamma \frac{\pi}{2}\right) \left[\left(c_\alpha \omega^\alpha\right)^2+\left(c_\beta \omega^\beta\right)^2 \right]+\left(c_\gamma \omega^\gamma\right)^2 \left[c_\alpha \omega^\alpha \sin\left(\alpha \frac{\pi}{2}\right)+c_\beta \omega^\beta \sin\left(\beta \frac{\pi}{2}\right) \right] + c_\alpha \omega^\alpha \cdot c_\beta\omega^\beta \cdot c_\gamma \omega^\gamma \left[\sin\left((\alpha-\beta-\gamma) \frac{\pi}{2}\right)+\sin\left((\beta-\alpha-\gamma) \frac{\pi}{2}\right) \right]}{\left(c_\alpha \omega^\alpha\right)^2+\left(c_\beta \omega^\beta\right)^2+\left(c_\gamma \omega^\gamma\right)^2+2c_\alpha \omega^\alpha \cdot c_\beta \omega^\beta \cos((\alpha-\beta)\frac{\pi}{2})+2c_\alpha \omega^\alpha \cdot c_\gamma \omega^\gamma \cos((\alpha-\gamma)\frac{\pi}{2})+2c_\beta \omega^\beta \cdot c_\gamma \omega^\gamma \cos((\beta-\gamma)\frac{\pi}{2})}
# ```

# # Fractional SLS (PT)

FractSLS_PT
#-

models = Vector{RheoModel}()

## plot moduli for varying α
for alpha in [0.1, 0.25, 0.5, 0.74, 0.9]
    
    push!(models, RheoModel(FractSLS_PT, (cₐ = 1, a = alpha, kᵦ = 1, kᵧ = 1)))
 
end

plotmodel(models)

# # Fractional Jeffreys (PT)

FractJeffreys_PT
#-

models = Vector{RheoModel}()

## plot moduli for varying β
for beta in [0.1, 0.25, 0.5, 0.74, 0.9]
    
    push!(models, RheoModel(FractJeffreys_PT, (ηₐ = 1, cᵦ = 1, β = beta, ηᵧ = 1)))
 
end

plotmodel(models, ymaxG = 0.5)

# # Standard Linear Solid (PT)

SLS_PT
#-

models = Vector{RheoModel}()

## plot moduli for varying kᵦ
for k in [1.0, 3.0, 5.0]
    
    push!(models, RheoModel(SLS_PT, (η = 1, kᵦ = k, kᵧ = 1)))
 
end

plotmodel(models)

# # Jeffreys (Zener)

Jeffreys_PT
#-

models = Vector{RheoModel}()

## plot moduli for varying ηₐ
for eta in [1.0, 5.0, 8.0]
    
    push!(models, RheoModel(Jeffreys_PT, (ηₐ = eta, k = 3, ηᵧ = 1)))
 
end

plotmodel(models)
