# # Springpot

using RHEOS

# By typing the name of the model, it is possible to visualise its graphical representation and its parameters.

Springpot

# ### Constitutive Equation
# ```math
# a = 1
# ```

# ### Relaxation Modulus
# ```math
# a = 1
# ```

# ### Creep Modulus
# ```math
# J(t)= frac{1}{c_{\beta} \Gamma (1 + \beta)} t^{\beta}
# ```

# ### Storage Modulus
# ```math
# a = 1
# ```

# ### Loss Modulus
# ```math
# a = 1
# ```

# # Spring

# When β = 0 the springpot specialises to a spring.

Spring

# # Dashpot

# When β = 1 the springpot specialises to a dashpot.

Dashpot

# ## Qualitative Behaviours of the Moduli

## include a helper function for plotting
include("assets/plothelper.jl") 

models = Vector{RheoModel}()

# Spring
push!(models, RheoModel(Spring, (k = 1.0,)))

# plot moduli for varying β
for beta in [0.2, 0.5, 0.8]
    
    push!(models, RheoModel(Springpot,(cᵦ = 1.0, β = beta)))
 
end

# Dashpot
push!(models, RheoModel(Dashpot, (η = 1.0,)))

plotmodel(models, ymaxG = 2.0)