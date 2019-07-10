function _nanremove()
    testarr = [1.0 1.0 NaN; NaN NaN 2.0; 3.0 3.0 3.0; NaN 4.0 4.0; 5.0 5.0 5.0]

    RHEOS.nanremove(testarr) == [3.0 3.0 3.0; 5.0 5.0 5.0]
end
@test _nanremove()

function _importdata_t_withnans()
    fildir = joinpath(@__DIR__, "testdata", "datawithnans.csv")

    rheodata1 = importdata(fildir; t_col = 1, σ_col = 2, ϵ_col = 3)
    test1 = (rheodata1.t == [5.0, 6.0] && rheodata1.σ == [10.0, 11.0] && rheodata1.ϵ == [15.0, 16.0])
    
    rheodata2 = importdata(fildir; t_col = 3, σ_col = 1, ϵ_col = 2)
    test2 = (rheodata2.t == [15.0, 16.0] && rheodata2.σ == [5.0, 6.0] && rheodata2.ϵ == [10.0, 11.0])

    rheodata3 = importdata(fildir; t_col = 2, σ_col = 3, ϵ_col = 1)
    test3 = (rheodata3.t == [10.0, 11.0] && rheodata3.σ == [15.0, 16.0] && rheodata3.ϵ == [5.0, 6.0])

    test1 && test2 && test3
end
@test _importdata_t_withnans()

function _importdata_ω_withnans()
    fildir = joinpath(@__DIR__, "testdata", "datawithnans.csv")

    rheodata1 = importdata(fildir; ω_col = 1, Gp_col = 2, Gpp_col = 3)
    test1 = (rheodata1.ω == [5.0, 6.0] && rheodata1.Gp == [10.0, 11.0] && rheodata1.Gpp == [15.0, 16.0])
    
    rheodata2 = importdata(fildir; ω_col = 3, Gp_col = 1, Gpp_col = 2)
    test2 = (rheodata2.ω == [15.0, 16.0] && rheodata2.Gp == [5.0, 6.0] && rheodata2.Gpp == [10.0, 11.0])

    rheodata3 = importdata(fildir; ω_col = 2, Gp_col = 3, Gpp_col = 1)
    test3 = (rheodata3.ω == [10.0, 11.0] && rheodata3.Gp == [15.0, 16.0] && rheodata3.Gpp == [5.0, 6.0])

    test1 #&& test2 && test3
end
@test _importdata_ω_withnans()