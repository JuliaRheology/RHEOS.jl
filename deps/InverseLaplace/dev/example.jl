using InverseLaplace

function test0()
    F = s -> 1/sqrt(s^2+1)
    t = [0.1, 1,10]
    N = 64
    sig =1.0
    b = 1.0
    f = weeks(F,t,N,sig,b)
    println(f)
    t1 = 0.1
    f = weeks(F,t1,N,sig,b)
    println(f)

    f = weekse(F,t1,N,sig,b)

end

function test1()
    F = s -> 1/sqrt(s^2+1)
    t = [0.1, 1,10]
    N = 64
    sig =1.0
    b = 1.0
    f = weeks(F,t,N,sig,b)
    fex = besselj0(t)
    println("abserr first: ",  abs(f-fex))
    (f, est) = weekse(F,t,N,sig,b)
    println("est error: ", est)
end

function test1a()
    F = s -> 1/sqrt(s^2+1)
    t = [0.1, 1,10]
    N = 64
    sig =1.0
    b = 1.0
    f = weeks(F,t,N,sig,b)
end


function test2()
    res = weeks( s -> s/(1+s^2) , [1.0,2.0,3.0],  18, 1.0,1.0)
    abs(res - [0.5403022953181169,-0.41614641727651525, -0.9899933071110818])
end

function test3()
    res = weekse( s -> s/(1+s^2) , [1.0,2.0,3.0],  18, 1.0,1.0)
    res
end

function test4(N=32)
    F(s) = exp(-4*sqrt(s))
    t = 1.0
    sig = 0.7
    b = 1.75
    res = weeks(F,[t],N,sig,b)[1]
    println(abs(res - 0.020675159969109663))
    fex = 2*exp(-4/t)/sqrt(pi*t^3)
    abserr = abs(res-fex)
    println("First abserr ", abserr)

    println("Difference from expected: ", abs(8.174615017609438e-6 - abserr))

    (so,bo) = wpar2(F,t,N,0,30,30)

    res = weeks(F,t,N,so,bo)

    abserr = abs(res-fex)
    println("Second abserr ", abserr)
    println("Difference from expected: ", abs(1.6406273631308643e-9  - abserr))    
end


function test4a(N=32)
    F(s) = exp(-4*sqrt(s))
    t = 1.0
    sig = 0.7
    b = 1.75
    local res
    @time (so,bo) = wpar2(F,t,N,0,30,30)
    @time res = weeks(F,t,N,so,bo)
    res
end

function timeweeks(F, t::AbstractVector, N, sig, b)
    t1 = @elapsed  a0 = real(wcoef(F,N,sig,b))
    a = a0[N+1:2*N]
    t2 = @elapsed L = laguer(a,2*b*t)
    t3 = @elapsed f = L .* exp((sig-b)*t)
    (t1,t2,t3)
end

function runweeks()
    F = s -> 1/sqrt(s^2+1)
    t = [0.1, 1,10]
    N = 80
    sig =1.0
    b = 1.0
    ts1 = 0.0
    ts2 = 0.0
    ts3 = 0.0
    for i in 1:3*10^2
        (t1,t2,t3) = timeweeks(F,t,N,sig,b)
        ts1 += t1
        ts2 += t2
        ts3 += t3
    end
#    (ts1,ts2,ts3)
    (1,ts2/ts1,ts3/ts1)    
end

function test5()
    F = s -> 1/sqrt(s^2+1)
    t = 100.0
    N = 100
    (so,bo) = wpar2(F,t,N,0,30,30)
    println("so $so, bo $bo")
    for t1 in 1//10:1//4:5
        println(t1)
        f = weeks(F,t1,N,so,bo)
        fex = besselj0(t1)
        ae  = abs(f-fex)
        ae2 = abs(ilt(F,BigFloat(t1), 64)- fex)
        @printf("t %.2f\t   ae %.4e    ae2 %.4e : fex %e ,  f %e \n", t1, ae, ae2, fex, f)
    end
end

