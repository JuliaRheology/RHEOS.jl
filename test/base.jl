# finite difference derivative
function deriv_test_float()

    x = collect(0.0:0.001:1.5)

    y = x.^2

    dy = 2*x

    dy_numeric = deriv(y, x)

    return sum((dy - dy_numeric).^2)

end

@test deriv_test_float() < 1e-5

function deriv_test_int_simple()

    x = collect(1:10)

    y = deriv(x)

    ysum = sum(y)

    return ysum

end

@test deriv_test_int_simple() == 10
@test typeof(deriv_test_int_simple()) <: Integer


