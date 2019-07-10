function _nanremove()
    testarr = [1.0 1.0 NaN; NaN NaN 2.0; 3.0 3.0 3.0; NaN 4.0 4.0; 5.0 5.0 5.0]

    RHEOS.nanremove(testarr) == [3.0 3.0 3.0; 5.0 5.0 5.0]
end
@test _nanremove()