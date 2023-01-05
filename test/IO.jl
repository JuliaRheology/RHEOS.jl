println("Testing IO.jl")
println("===============================================")



function _nanremove_withnans()
    testarr = [1.0 1.0 NaN; NaN NaN 2.0; 3.0 3.0 3.0; NaN 4.0 4.0; 5.0 5.0 5.0]

    RHEOS.nanremove(testarr) == [3.0 3.0 3.0; 5.0 5.0 5.0]
end
@test _nanremove_withnans()

function _nanremove_nonans()
    testarr = [1.0 1.0 1.0; 2.0 2.0 2.0; 3.0 3.0 3.0; 4.0 4.0 4.0; 5.0 5.0 5.0]

    RHEOS.nanremove(testarr) == [1.0 1.0 1.0; 2.0 2.0 2.0; 3.0 3.0 3.0; 4.0 4.0 4.0; 5.0 5.0 5.0]
end
@test _nanremove_nonans()

function _importcsv_t_withnans()
    fildir = joinpath(@__DIR__, "testdata", "datawithnans.csv")

    rheodata1 = importcsv(fildir; t_col = 1, σ_col = 2, ϵ_col = 3)
    test1 = (rheodata1.t == [5.0, 6.0] && rheodata1.σ == [10.0, 11.0] && rheodata1.ϵ == [15.0, 16.0])

    rheodata2 = importcsv(fildir; t_col = 3, σ_col = 1, ϵ_col = 2)
    test2 = (rheodata2.t == [15.0, 16.0] && rheodata2.σ == [5.0, 6.0] && rheodata2.ϵ == [10.0, 11.0])

    rheodata3 = importcsv(fildir; t_col = 2, σ_col = 3, ϵ_col = 1)
    test3 = (rheodata3.t == [10.0, 11.0] && rheodata3.σ == [15.0, 16.0] && rheodata3.ϵ == [5.0, 6.0])

    test1 && test2 && test3
end
@test _importcsv_t_withnans()

function _importcsv_ω_withnans()
    fildir = joinpath(@__DIR__, "testdata", "datawithnans.csv")

    rheodata1 = importcsv(fildir; ω_col = 1, Gp_col = 2, Gpp_col = 3)
    test1 = (rheodata1.ω == [5.0, 6.0] && rheodata1.Gp == [10.0, 11.0] && rheodata1.Gpp == [15.0, 16.0])

    rheodata2 = importcsv(fildir; ω_col = 3, Gp_col = 1, Gpp_col = 2)
    test2 = (rheodata2.ω == [15.0, 16.0] && rheodata2.Gp == [5.0, 6.0] && rheodata2.Gpp == [10.0, 11.0])

    rheodata3 = importcsv(fildir; ω_col = 2, Gp_col = 3, Gpp_col = 1)
    test3 = (rheodata3.ω == [10.0, 11.0] && rheodata3.Gp == [15.0, 16.0] && rheodata3.Gpp == [5.0, 6.0])

    test1 && test2 && test3
end
@test _importcsv_ω_withnans()

function _importcsv_t_partial()
    fildir = joinpath(@__DIR__, "testdata", "datapartial.csv")

    rheodata1 = importcsv(fildir; t_col = 1)
    test1 = (rheodata1.t == [1.0, 2.0, 3.0, 4.0, 5.0])

    rheodata2 = importcsv(fildir; t_col = 1, σ_col = 2)
    test2 = (rheodata2.t == [1.0, 2.0, 3.0, 4.0, 5.0] && rheodata2.σ == [10.0, 20.0, 30.0, 40.0, 50.0])

    rheodata3 = importcsv(fildir; t_col = 1, ϵ_col = 1)
    test3 = (rheodata3.t == [1.0, 2.0, 3.0, 4.0, 5.0] && rheodata3.ϵ == [10.0, 20.0, 30.0, 40.0, 50.0])

    test1 && test2 && test3
end
@test _importcsv_t_withnans()

function _exportcsv_timestress()
    fildir = joinpath(@__DIR__, "testdata", "datapartial.csv")

    rheotimedataIN = importcsv(fildir; t_col = 1, σ_col = 2)

    # default column ordering for full time data is (t, σ)
    # testdir = joinpath(@__DIR__, "testdata", "testloop.csv")
    testdir = tempname()
    exportcsv(rheotimedataIN, testdir)
    rheotimedataOUT = importcsv(testdir; t_col = 1, σ_col = 2)
    rm(testdir)
    rheotimedataIN.σ == rheotimedataOUT.σ && rheotimedataIN.ϵ == rheotimedataOUT.ϵ && rheotimedataIN.t == rheotimedataOUT.t
end
@test _exportcsv_timestress()

function _exportcsv_timestrain()
    fildir = joinpath(@__DIR__, "testdata", "datapartial.csv")

    rheotimedataIN = importcsv(fildir; t_col = 1, ϵ_col = 2)

    # default column ordering for full time data is (t, ϵ)
    testdir = tempname()
    exportcsv(rheotimedataIN, testdir)
    rheotimedataOUT = importcsv(testdir; t_col = 1, ϵ_col = 2)
    rm(testdir)
    rheotimedataIN.σ == rheotimedataOUT.σ && rheotimedataIN.ϵ == rheotimedataOUT.ϵ && rheotimedataIN.t == rheotimedataOUT.t
end
@test _exportcsv_timestrain()

function _exportcsv_timefull()
    fildir = joinpath(@__DIR__, "testdata", "datanonans.csv")

    rheotimedataIN = importcsv(fildir; t_col = 1, σ_col = 2, ϵ_col = 3)

    # default column ordering for full time data is (t, σ, ϵ)
    testdir = tempname()
    exportcsv(rheotimedataIN, testdir)
    rheotimedataOUT = importcsv(testdir; t_col = 1, ϵ_col = 3, σ_col = 2)
    rm(testdir)
    rheotimedataIN.σ == rheotimedataOUT.σ && rheotimedataIN.ϵ == rheotimedataOUT.ϵ && rheotimedataIN.t == rheotimedataOUT.t
end
@test _exportcsv_timefull()

function _exportcsv_freqfull()
    fildir = joinpath(@__DIR__, "testdata", "datanonans.csv")

    rheofreqdataIN = importcsv(fildir; ω_col = 1, Gp_col = 2, Gpp_col = 3)

    # default column ordering for full time data is (ω, Gp, Gpp)
    testdir = tempname()
    exportcsv(rheofreqdataIN, testdir)
    rheofreqdataOUT = importcsv(testdir; ω_col = 1, Gp_col = 2, Gpp_col = 3)
    rm(testdir)
    rheofreqdataIN.ω == rheofreqdataOUT.ω && rheofreqdataIN.Gp == rheofreqdataOUT.Gp && rheofreqdataIN.Gpp == rheofreqdataOUT.Gpp
end
@test _exportcsv_freqfull()

# function _saveload_timefull()
#     fildir = joinpath(@__DIR__, "testdata", "datanonans.csv")

#     rheotimedataIN = importcsv(fildir; t_col = 1, σ_col = 2, ϵ_col = 3)

#     # default column ordering for full time data is (σ, ϵ, t)
    # testdir = tempname()
#     savedata(rheotimedataIN, testdir)
#     rheotimedataOUT = loaddata(testdir)

    # rm(testdir)
#     rheotimedataIN.σ == rheotimedataOUT.σ && rheotimedataIN.ϵ == rheotimedataOUT.ϵ && rheotimedataIN.t == rheotimedataOUT.t
# end
# @test _saveload_timefull()

# function _saveload_freqfull()
#     fildir = joinpath(@__DIR__, "testdata", "datanonans.csv")

#     rheofreqdataIN = importcsv(fildir; ω_col = 1, Gp_col = 2, Gpp_col = 3)

#     # default column ordering for full time data is (Gp, Gpp, ω)
    # testdir = tempname()
#     savedata(rheofreqdataIN, testdir)
#     rheofreqdataOUT = loaddata(testdir)

    # rm(testdir)
#     rheofreqdataIN.ω == rheofreqdataOUT.ω && rheofreqdataIN.Gp == rheofreqdataOUT.Gp && rheofreqdataIN.Gpp == rheofreqdataOUT.Gpp
# end
# @test _saveload_freqfull()
