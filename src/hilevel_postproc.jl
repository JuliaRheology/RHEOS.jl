#!/usr/bin/env julia

"""
    fiteval(self::RheologyData, modelname::String; singularity = false)

Show plot of data vs. fitted data for specified model.
"""
function fiteval(self::RheologyData, modelname::String; singularity = false)

    # params
    params = self.fittedmodels[modelname]

    # modulus function
    model = moduli(modelname, self.test_type)

    # get data
    if self.test_type == "strlx"
        measured = self.σ
        prescribed_dot = self.dϵ
    elseif self.test_type == "creep"
        measured = self.ϵ
        prescribed_dot = self.dσ
    end

    # get fit
    fitted = boltzconvolve(model, self.t, deriv(self.t), params, prescribed_dot; singularity = singularity)

    # print params
    println(modelname, " fit: ", self.fittedmodels[modelname])

    # stress subplot
    if singularity
        plot(self.t[1:end], measured)
        plot(self.t[2:end], fitted, "--")
        show()
    else
        plot(self.t, measured)
        plot(self.t, fitted, "--")
        show()
    end

end

# function that saves data to file

# add legends to plot
