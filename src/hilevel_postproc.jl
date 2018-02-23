#!/usr/bin/env julia

"""
    fiteval(self::RheologyData, modelname::String, test_type::String; singularity = false)

Show plot of data vs. fitted data for specified model.
"""
function fiteval(self::RheologyData, modelname::String, test_type::String; singularity = false)

    # params
    params = self.fittedmodels[modelname]

    # modulus function
    model = moduli(modelname, test_type)

    # get data
    if test_type == "strlx"
        measured = self.σᵦ
        prescribed_dot = self.dϵᵦ
    elseif test_type == "creep"
        measured = self.ϵᵦ
        prescribed_dot = self.dσᵦ
    end

    # get fit
    fitted = boltzconvolve(model, self.tᵦ, deriv(self.tᵦ), params, prescribed_dot; singularity = singularity)

    # print params
    println(modelname, " fit: ", self.fittedmodels[modelname])

    # stress subplot
    if singularity
        plot(self.tᵦ[1:end], measured)
        plot(self.tᵦ[2:end], fitted, "--")
        show()
    else
        plot(self.tᵦ, measured)
        plot(self.tᵦ, fitted, "--")
        show()
    end

end

# function that saves data to file

# add legends to plot
