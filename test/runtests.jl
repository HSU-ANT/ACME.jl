# Copyright 2015, 2016, 2017, 2018, 2019, 2020, 2021 Martin Holters
# See accompanying license file.

include("checklic.jl")

using ACME
using Compat: evalpoly, only
using FFTW: rfft
using ProgressMeter: Progress, next!
using SparseArrays: sparse, spzeros
using Test: @test, @test_broken, @test_logs, @test_throws, @testset

@testset "topomat" begin
    tv, ti = ACME.topomat(sparse([1 -1 1; -1 1 -1]))
    @test tv*ti'==spzeros(2,1)

    # Pathological cases for topomat:
    # two nodes, one loop branch (short-circuited) -> voltage==0, current arbitrary
    @test ACME.topomat(spzeros(Int, 2, 1)) == (hcat([1]), spzeros(0, 1))
    # two nodes, one branch between them -> voltage arbitrary, current==0
    @test ACME.topomat(sparse([1,2], [1,1], [1,-1])) == (spzeros(0, 1), hcat([1]))
end

@testset "LinearSolver" begin
    solver = ACME.LinearSolver(3)
    A = [1.0 0.5 0.4; 2.0 4.0 1.7; 4.0 7.0 9.1]
    @test ACME.setlhs!(solver, A)
    x = rand(3)
    y = similar(x)
    ACME.solve!(solver, y, x)
    @test A*y ≈ x
    y = copy(x)
    ACME.solve!(solver, y, y)
    @test A*y ≈ x
    @test_throws DimensionMismatch ACME.setlhs!(solver, zeros(2, 3))
    @test_throws DimensionMismatch ACME.setlhs!(solver, zeros(3, 4))
    @test_throws DimensionMismatch ACME.setlhs!(solver, zeros(4, 4))
    @test_throws DimensionMismatch ACME.solve!(solver, zeros(2), zeros(3))
    @test_throws DimensionMismatch ACME.solve!(solver, zeros(3), zeros(4))
    @test_throws DimensionMismatch ACME.solve!(solver, zeros(4), zeros(4))
    @test !ACME.setlhs!(solver, zeros(3,3))
end

@testset "simple circuits" begin
    @testset "empty circuit" begin
        circ = @circuit begin end
        model = DiscreteModel(circ, 1)
        @test run!(model, zeros(0, 20)) == zeros(0, 20)
    end

    @testset "only one resistor" begin
        circ = @circuit begin
            r = resistor(0), [1] ⟷ [2]
        end
        model = DiscreteModel(circ, 1)
        @test run!(model, zeros(0, 20)) == zeros(0, 20)
    end

    # simple circuit: resistor and diode in series, driven by constant voltage,
    # chosen such that a prescribe current flows
    @testset "resistor-diode" begin
        i = 1e-3
        r=10e3
        is=1e-12
        v_r = i*r
        v_d = 25e-3 * log(i/is+1)
        circ = @circuit begin
            vsrc = voltagesource(v_r + v_d), [+] ⟷ "supply voltage",[-] ⟷ gnd
            r1 = resistor(r)
            d = diode(is=is), [-] ⟷ gnd, [+] ⟷ r1[2]
            vprobe = voltageprobe(), [-] ⟷ gnd, [+] ⟷ r1[2]
            r1[1] ⟷ "supply voltage"
        end
        model = DiscreteModel(circ, 1)
        y = run!(model, zeros(0, 1))
        @test y[1] ≈ v_d
    end
end

@testset "circuit with reused resdef" begin
    @test_logs (:warn, "redefinition of `r1`") macroexpand(
        @__MODULE__,
        :(
            @circuit begin
                r1 = resistor(100), [1] ↔ gnd
                r2 = resistor(100), [1] ↔ r1[2]
                r1 = resistor(100), [1] ↔ r2[2]
            end
        )
    )
end

@testset "circuit manipulation" begin
    @testset "programmatic reconnection" begin
        circ = @circuit begin
            r1 = resistor(10)
            r2 = resistor(100), [1] ⟷ r1[1], [2] ⟷ r1[2]
            src = voltagesource(1), [-] ⟷ r1[2]
            probe = currentprobe(), [+] ⟷ src[+],
                                    [-] ⟷ r1[1]
        end
        model = DiscreteModel(circ, 1)
        @test run!(model, zeros(0, 1))[1,1] ≈ 1/10 + 1/100
        disconnect!(circ, (:r2, 1))
        model = DiscreteModel(circ, 1)
        @test run!(model, zeros(0, 1))[1,1] ≈ 1/10
        disconnect!(circ, (:r1, 2))
        model = DiscreteModel(circ, 1)
        @test run!(model, zeros(0, 1))[1,1] ≈ 0
        connect!(circ, (:r1, 2), (:r2, 1))
        model = DiscreteModel(circ, 1)
        @test run!(model, zeros(0, 1))[1,1] ≈ 1/(10+100)
    end

    @testset "element deletion" begin
        circ = @circuit begin
            r1 = resistor(10)
        end
        r2_des = add!(circ, resistor(100))
        add!(circ, :r3, resistor(470))
        r4_des = add!(circ, resistor(1000))
        add!(circ, :src, voltagesource(1))
        add!(circ, :probe, currentprobe())
        connect!(circ, (:src, :+), (:probe, :+))
        connect!(circ, (:probe, :-), (:r1, 1), (r2_des, 1), (:r3, 1), (r4_des, 1))
        connect!(circ, (:src, :-), (:r1, 2), (r2_des, 2), (:r3, 2), (r4_des, 2))
        model = DiscreteModel(circ, 1)
        @test run!(model, zeros(0, 1))[1,1] ≈ 1/10 + 1/100 + 1/470 + 1/1000
        delete!(circ, :r1)
        model = DiscreteModel(circ, 1)
        @test run!(model, zeros(0, 1))[1,1] ≈ 1/100 + 1/470 + 1/1000
        delete!(circ, r4_des)
        model = DiscreteModel(circ, 1)
        @test run!(model, zeros(0, 1))[1,1] ≈ 1/100 + 1/470
        delete!(circ, :r3)
        model = DiscreteModel(circ, 1)
        @test run!(model, zeros(0, 1))[1,1] ≈ 1/100
        delete!(circ, r2_des)
        model = DiscreteModel(circ, 1)
        @test run!(model, zeros(0, 1))[1,1] ≈ 0
    end
end

@testset "faulty circuits" begin
    @testset "indeterminate output" begin
        circ = @circuit begin
            r = resistor(0)
            probe = currentprobe(), [+] ⟷ r[1], [-] ⟷ r[2]
        end
        @test_logs (:warn, "Model output depends on indeterminate quantity") DiscreteModel(circ, 1)
    end

    @testset "no solution" begin
        circ = @circuit begin
            d = diode()
            src = currentsource(), [+] ↔ d[+], [-] ↔ d[-]
            probe = voltageprobe(), [+] == d[+], [-] == d[-]
        end
        model = DiscreteModel(circ, 1)
        @test ACME.nn(model) == 1
        y = run!(model, [1.0 1.0])
        @test size(y) == (1, 2)
        @test y[1,1] == y[1,2]
        @test_throws ErrorException run!(model, hcat([Inf]))
        @test(size(@test_logs((:warn, "Failed to converge while solving non-linear equation."), run!(model, hcat([-1.0])))) == (1, 1))
    end
end

@testset "KDTree" begin
    @testset "size = 4×$num" for num = 1:50
        let ps = rand(4, num)
            t = ACME.KDTree(ps)
            for i in 1:size(ps)[2]
                idx = ACME.indnearest(t, ps[:,i])
                @test ps[:,i] == ps[:,idx]
            end
        end
    end

    @testset "size = 5×10000" begin
        ps = rand(6, 10000)
        t = ACME.KDTree(ps)
        p = rand(6)
        best_p = ps[:,argmin(vec(sum(abs2, ps .- p, dims=1)))]
        idx = ACME.indnearest(t, p)
        @test sum(abs2, p - best_p) ≈ sum(abs2, p - ps[:, idx])
    end
end

@testset "HomotopySolver" begin
    nleq = ACME.ParametricNonLinEq((res, J, scratch, z) ->
    let p=scratch[1], Jp=scratch[2]
        res[1] = z[1]^2 - 1 + p[1]
        J[1,1] = 2*z[1]
        Jp[1,1] = 1
    end, 1, 1)
    solver = ACME.HomotopySolver{ACME.SimpleSolver}(nleq, [0.0], [1.0])
    ACME.solve(solver, [-0.5 + rand()])
    @test ACME.hasconverged(solver)
    ACME.solve(solver, [1.5 + rand()])
    @test !ACME.hasconverged(solver)
end

@testset "gensolve/rank_factorize" begin
    a = Rational{BigInt}[1 1 1; 1 1 2; 1 2 1; 1 2 2; 2 1 1; 2 1 2]
    b = Rational{BigInt}[1 2 3 4 5 6; 6 5 4 3 2 1; 1 0 1 0 1 0]
    nullspace = ACME.gensolve(sparse(a'), spzeros(Rational{BigInt}, size(a, 2), 0))[2]
    @test nullspace'*a == spzeros(3, 3)
    c, f = ACME.rank_factorize(sparse(a * b))
    @test c*f == a*b
end

@testset "redcue_pdims" begin
    a = Rational{BigInt}[-1 -1 -4 -3 0 -1; 2 -1 -5 3 -4 0; -2 2 -5 -2 5 1;
                         -5 4 -3 0 5 -5; 4 3 0 -1 0 2; 0 -3 -4 -4 -3 4]
    b = hcat(Rational{BigInt}[1; 2; 3; -2; -1; 0])
    c = Rational{BigInt}[4 2 -1; -1 -3 0; -3 5 3; 0 0 0; -4 -1 -1; -1 -1 5]
    dy = Rational{BigInt}[1 2 3 -2 -1 0]
    ey = hcat(Rational{BigInt}[5])
    fy = Rational{BigInt}[-2 -1 3]
    p = Rational{BigInt}[1 1 1; 1 1 2; 1 2 1; 1 2 2; 2 1 1; 2 1 2]
    dq = Rational{BigInt}[1 2 3 4 5 6; 6 5 4 3 2 1; 1 0 1 0 1 0]
    eq = hcat(Rational{BigInt}[1; 2; 3])
    fq = Rational{BigInt}[1 0 0; 10 0 0; 0 1 0; 0 10 0; 0 0 1; 0 0 10]

    for zxin in (zeros(Rational{BigInt}, 3, 6), Rational{BigInt}[1 2 0 0 2 1; 0 1 2 2 0 1; 0 0 1 0 1 1]),
        zuin in (zeros(Rational{BigInt}, 3, 1), hcat(Rational{BigInt}[1; 2; -1]))
        dq_full = p * dq + fq * zxin
        eq_full = p * eq + fq * zuin
        mats = Dict{Symbol,Array}(:a => a, :b => b, :c => c, :dy => dy, :ey => ey,
            :fy => fy, :dq_full => dq_full, :eq_full => eq_full, :fq => fq)
        mats[:dq_fulls]=Matrix[mats[:dq_full]]
        mats[:eq_fulls]=Matrix[mats[:eq_full]]
        mats[:fqprev_fulls]=Matrix[mats[:eq_full]]
        mats[:fqs]=Matrix[mats[:fq]]
        ACME.reduce_pdims!(mats)
        @test size(mats[:pexps][1], 2) == 3
        @test mats[:pexps][1] * mats[:dqs][1] == mats[:dq_fulls][1]
        @test mats[:pexps][1] * mats[:eqs][1] == mats[:eq_fulls][1]
        zx = (fq'fq) \ fq' * (dq_full - mats[:dq_fulls][1])
        zu = (fq'fq) \ fq' * (eq_full - mats[:eq_fulls][1])
        @test mats[:a] == a - c*zx
        @test mats[:b] == b - c*zu
        @test mats[:dy] == dy - fy*zx
        @test mats[:ey] == ey - fy*zu
    end

    # TODO: Also check dq, eq, fqprev for a test case with desomposed non-linearity
end

@testset "nonlinearity decomposition" begin
    circ = @circuit begin
        src1 = voltagesource()
        probe1 = currentprobe()
        d1 = diode(), [+] ⟷ src1[+]
        d2 = diode(), [+] ⟷ d1[-], [-] ⟷ probe1[+]
        probe1[-] ⟷ src1[-]
    end
    add!(circ, :src2, voltagesource())
    add!(circ, :probe2, currentprobe())
    add!(circ, :d3, diode())
    connect!(circ, (:src2, :+), (:d3, :+))
    connect!(circ, (:d3, :-), (:probe2, :+))
    connect!(circ, (:probe2, :-), (:src2, :-))
    model = DiscreteModel(circ, 1, decompose_nonlinearity=false)
    y = run!(model, hcat([2.0; 1.0]))
    @test ACME.nn(model, 1) == 3
    @test y[1] ≈ 1e-12*(exp(1/25e-3)-1)
    @test y[2] ≈ 1e-12*(exp(1/25e-3)-1)
    model = DiscreteModel(circ, 1)
    y = run!(model, hcat([2.0; 1.0]))
    # single diode is extracted first, although it was added last
    @test ACME.nn(model, 1) == 1
    @test ACME.nn(model, 2) == 2
    @test y[1] ≈ y[2] ≈ 1e-12*(exp(1/25e-3)-1)
end

@testset "composite_element" begin
    subcirc1() = composite_element(@circuit(begin
       r1 = resistor(100e3)
       r2 = resistor(1e3)
       r1[2] == r2[1]
       src = voltagesource(5), [+] == r1[1], [-] == r2[2]
   end), pinmap=Dict(1 => (:r2, 1), 2 => (:r2, 2)))
    circ = @circuit begin
        U = subcirc1()
        J = voltageprobe(gp=2), [+] == U[1], [-] == U[2]
    end
    model = DiscreteModel(circ, 1//44100)
    y = run!(model, zeros(0, 100))
    refcirc = @circuit begin
        r1 = resistor(100e3)
        r2 = resistor(1e3)
        r1[2] == r2[1]
        src = voltagesource(5), [+] == r1[1], [-] == r2[2]
        J = voltageprobe(gp=2), [+] == r2[1], [-] == r2[2]
    end
    refmodel = DiscreteModel(refcirc, 1//44100)
    yref = run!(refmodel, zeros(0, 100))
    @test y ≈ yref

    subcirc2() = composite_element(@circuit(begin
       r1 = resistor(100e3)
       r2 = resistor(1e3)
       r1[2] == r2[1]
       src = voltagesource(), [+] == r1[1], [-] == r2[2]
   end), pinmap=Dict(1 => (:r2, 1), 2 => (:r2, 2)))
    circ = @circuit begin
        U = subcirc2()
        J = voltageprobe(gp=2), [+] == U[1], [-] == U[2]
    end
    model = DiscreteModel(circ, 1//44100)
    y = run!(model, 5*ones(1, 100))
    yref = run!(refmodel, zeros(0, 100))
    @test y ≈ yref

    subcirc3() = composite_element(@circuit(begin
       r1 = resistor(100e3)
       r2 = resistor(1e3)
       c = capacitor(1e-6), [1] == r2[1], [2] == r2[2]
       r1[2] == r2[1]
       src = voltagesource(5), [+] == r1[1], [-] == r2[2]
   end), pinmap=Dict(1 => (:r2, 1), 2 => (:r2, 2)))
    circ = @circuit begin
        U = subcirc3()
        J = voltageprobe(gp=2), [+] == U[1], [-] == U[2]
    end
    model = DiscreteModel(circ, 1//44100)
    y = run!(model, zeros(0, 100))
    refcirc = @circuit begin
        r1 = resistor(100e3)
        r2 = resistor(1e3)
        c = capacitor(1e-6), [1] == r2[1], [2] == r2[2]
        r1[2] == r2[1]
        src = voltagesource(5), [+] == r1[1], [-] == r2[2]
        J = voltageprobe(gp=2), [+] == r2[1], [-] == r2[2]
    end
    refmodel = DiscreteModel(refcirc, 1//44100)
    yref = run!(refmodel, zeros(0, 100))
    @test y ≈ yref

    subcirc4() = composite_element(@circuit(begin
       r1 = resistor(100e3)
       r2 = resistor(1e3)
       c = capacitor(1e-6), [1] == r2[1], [2] == r2[2]
       d = diode(), [+] == r2[1], [-] == r2[2]
       r1[2] == r2[1]
       src = voltagesource(5), [+] == r1[1], [-] == r2[2]
   end), pinmap=Dict(1 => (:r2, 1), 2 => (:r2, 2)))
    circ = @circuit begin
        U = subcirc4()
        J = voltageprobe(gp=2), [+] == U[1], [-] == U[2]
    end
    model = DiscreteModel(circ, 1//44100)
    y = run!(model, zeros(0, 100))
    refcirc = @circuit begin
        r1 = resistor(100e3)
        r2 = resistor(1e3)
        c = capacitor(1e-6), [1] == r2[1], [2] == r2[2]
        d = diode(), [+] == r2[1], [-] == r2[2]
        r1[2] == r2[1]
        src = voltagesource(5), [+] == r1[1], [-] == r2[2]
        J = voltageprobe(gp=2), [+] == r2[1], [-] == r2[2]
    end
    refmodel = DiscreteModel(refcirc, 1//44100)
    yref = run!(refmodel, zeros(0, 100))
    @test y ≈ yref
end

@testset "sources/probes" begin
# sources and probes with internal resistance/conductance
    circ = @circuit begin
        src = currentsource(100e-3, gp=1//100000)
        probe = voltageprobe(), [+] ⟷ src[+], [-] ⟷ src[-]
    end
    model = DiscreteModel(circ, 1)
    @test run!(model, zeros(0,1)) ≈ [100000*100e-3]

    circ = @circuit begin
        src = currentsource(gp=1//100000)
        probe = voltageprobe(), [+] ⟷ src[+], [-] ⟷ src[-]
    end
    model = DiscreteModel(circ, 1)
    @test run!(model, hcat([100e-3])) ≈ [100000*100e-3]

    circ = @circuit begin
        src = currentsource(100e-3)
        probe = voltageprobe(gp=1//100000), [+] ⟷ src[+], [-] ⟷ src[-]
    end
    model = DiscreteModel(circ, 1)
    @test run!(model, zeros(0,1)) ≈ [100000*100e-3]

    circ = @circuit begin
        src = voltagesource(10, rs=100000)
        probe = currentprobe(), [+] ⟷ src[+], [-] ⟷ src[-]
    end
    model = DiscreteModel(circ, 1)
    @test run!(model, zeros(0,1)) ≈ [10/100000]

    circ = @circuit begin
        src = voltagesource(rs=100000)
        probe = currentprobe(), [+] ⟷ src[+], [-] ⟷ src[-]
    end
    model = DiscreteModel(circ, 1)
    @test run!(model, hcat([10.0])) ≈ [10/100000]

    circ = @circuit begin
        src = voltagesource(10)
        probe = currentprobe(rs=100000), [+] ⟷ src[+], [-] ⟷ src[-]
    end
    model = DiscreteModel(circ, 1)
    @test run!(model, zeros(0,1)) ≈ [10/100000]
end

@testset "Jiles-Atherton inductor/transformer" begin
    @testset "Jiles-Atherton inductor" begin
        circ = @circuit begin
            Jin = voltagesource()
            Jout1 = currentprobe(), [+] ↔ Jin[+]
            Jout2 = currentprobe(), [+] ↔ Jin[+]
            L_JA = inductor(Val{:JA}), [1] ↔ Jout1[-], [2] ↔ Jin[-]
            # default parameters approximately linearize as 174mH
            L_lin = inductor(174e-3), [1] ↔ Jout2[-], [2] ↔ Jin[-]
        end
        model = DiscreteModel(circ, 1//44100)
        # starting non-magnetized, the inductor first behaves sub-linear
        y = run!(model, fill(0.1, 1, 750))
        @test isapprox(y[1,1:9], y[2,1:9], rtol=1e-2) # almost linear at first
        @test all(y[1,:] .< y[2,:])
        # until it reaches saturation, where it becomes super-linear
        run!(model, fill(0.1, 1, 500))
        y = run!(model, fill(0.1, 1, 750))
        @test all(y[1,:] .> y[2,:])
        # due to hysteresis, applying the negative voltage for the same time
        # drives current below zero
        y = run!(model, fill(-0.1, 1, 2000))
        @test y[1,end] < -2e-3
        # with voltage at zero (i.e. shorted), the current stays constant
        y = run!(model, zeros(1, 1000))
        @test y[1,1] < -2e-3 && all(y[:,1] .≈ y)
    end
    @testset "Jiles-Atherton transformer" begin
        circ = @circuit begin
            Jin = voltagesource()
            R1 = resistor(10), [1] ↔ Jin[+]
            R2 = resistor(10), [1] ↔ Jin[+]
            T_JA = transformer(Val{:JA}; ns = [10, 100]), [1] ↔ R1[2], [2] ↔ Jin[-]
            # approximate small-scale equivalent
            T_lin = transformer(330e-6, 33e-3), [primary1] ↔ R2[2], [primary2] ↔ Jin[-]
            Jout1 = voltageprobe(;gp=1e-3), [+] ↔ T_JA[3], [-] ↔ T_JA[4]
            Jout2 = voltageprobe(;gp=1e-3), [+] ↔ T_lin[secondary1], [-] ↔ T_lin[secondary2]
        end
        model = DiscreteModel(circ, 1//44100)
        u = [sin(2pi*1000/44100*n) for n in (0:499)]'
        # almost linear for small input
        y = run!(model, 0.001*u)[:,200:end]
        @test isapprox(y[1,:], y[2,:]; rtol = 1.2e-3)
        y = run!(model, 0.002*u)[:,200:end]
        @test isapprox(y[1,:], y[2,:]; rtol = 1.2e-3)
        # not at all linear for large input
        y = run!(model, 10*u)[200:end]
        @test !isapprox(y[1,:], y[2,:]; rtol = 0.75)
    end
end

@testset "BJT" begin
    isc=1e-6
    ise=2e-6
    ηc=1.1
    ηe=1.0
    βf=100
    βr=10
    @testset "Ebers-Moll $typ" for (typ, ib) in ((:npn, 1e-3), (:pnp, -1e-3))
        circ = @circuit begin
            t = bjt(typ, isc=isc, ise=ise, ηc=ηc, ηe=ηe, βf=βf, βr=βr)
            isrc = currentsource(), [+] ⟷ t[base]
            vsrc = voltagesource(), [-] ⟷ isrc[-]
            veprobe = voltageprobe(), [+] ⟷ t[base], [-] ⟷ isrc[-]
            vcprobe = voltageprobe(), [+] ⟷ t[base], [-] ⟷ vsrc[+]
            ieprobe = currentprobe(), [+] ⟷ t[emitter], [-] ⟷ isrc[-]
            icprobe = currentprobe(), [+] ⟷ t[collector], [-] ⟷ vsrc[+]
        end
        model = DiscreteModel(circ, 1)
        N = 100
        output = run!(model, [range(0, stop=ib, length=N)'; range(1, stop=-1, length=N÷2)' range(-1, stop=1, length=N÷2)'])
        if typ == :pnp
            output = -output
        end
        for n in 1:N
            ve, vc, ie, ic = output[:,n]
            @test isapprox(ie, ise*(exp(ve/(ηe*25e-3))-1) - βr/(1+βr)*isc*(exp(vc/(ηc*25e-3))-1), atol=1e-10)
            @test isapprox(ic, -βf/(1+βf)*ise*(exp(ve/(ηe*25e-3))-1) + isc*(exp(vc/(ηc*25e-3))-1), atol=1e-10)
        end
    end
    ηcl=1.2
    ηel=1.3
    prog = Progress(2^9, "Testing Gummel-Poon model: ")
    @testset "Gummel-Poon $typ" for ile in (0, 50e-9), ilc in (0, 100e-9),
            ηcl in (ηc, 1.2), ηel in (ηe, 1.1),
            vaf in (Inf, 10), var in (Inf, 50),
            ikf in (Inf, 50e-3), ikr in (Inf, 500e-3),
            (typ, ib) in ((:npn, 1e-3), (:pnp, -1e-3))
        circ = @circuit begin
            t = bjt(typ, isc=isc, ise=ise, ηc=ηc, ηe=ηe, βf=βf, βr=βr, ile=ile,
                    ilc=ilc, ηcl=ηcl, ηel=ηel, vaf=vaf, var=var, ikf=ikf, ikr=ikr)
            isrc = currentsource(), [+] ⟷ t[base]
            vsrc = voltagesource(), [-] ⟷ isrc[-]
            veprobe = voltageprobe(), [+] ⟷ t[base], [-] ⟷ isrc[-]
            vcprobe = voltageprobe(), [+] ⟷ t[base], [-] ⟷ vsrc[+]
            ieprobe = currentprobe(), [+] ⟷ t[emitter], [-] ⟷ isrc[-]
            icprobe = currentprobe(), [+] ⟷ t[collector], [-] ⟷ vsrc[+]
        end
        model = DiscreteModel(circ, 1)
        N = 100
        output = run!(model, [range(0, stop=ib, length=N)'; range(1, stop=-1, length=N÷2)' range(-1, stop=1, length=N÷2)'])
        if typ == :pnp
            output = -output
        end
        for n in 1:N
            ve, vc, ie, ic = output[:,n]
            i_f = βf/(1+βf)*ise*(exp(ve/(ηe*25e-3))-1)
            i_r = βr/(1+βr)*isc*(exp(vc/(ηc*25e-3))-1)
            icc = (2*(1-ve/var-vc/vaf))/(1+sqrt(1+4(i_f/ikf+i_r/ikr))) * (i_f - i_r)
            ibe = 1/βf*i_f + ile*(exp(ve/(ηel*25e-3))-1)
            ibc = 1/βr*i_r + ilc*(exp(vc/(ηcl*25e-3))-1)
            @test isapprox(ie, icc + ibe, atol=1e-10)
            @test isapprox(ic, -icc + ibc, atol=1e-10)
        end
        next!(prog)
    end

    rb=100
    re=10
    rc=20
    @testset "internal resistances ($typ)" for (typ, ib, vce) in ((:npn, 1e-3, 1), (:pnp, -1e-3, -1))
        circ = @circuit begin
            t1 = bjt(typ)
            rbref = resistor(rb)
            rcref = resistor(rc)
            reref = resistor(re)
            isrc1 = currentsource(ib)
            vscr1 = voltagesource(vce)
            veprobe1 = voltageprobe()
            vcprobe1 = voltageprobe()
            ieprobe1 = currentprobe()
            icprobe1 = currentprobe()
            t1[base] ⟷ rbref[1]
            rbref[2] ⟷ isrc1[+] ⟷ veprobe1[+] ⟷ vcprobe1[+]
            t1[collector] ⟷ rcref[1]
            rcref[2] ⟷ icprobe1[+]
            vcprobe1[-] ⟷ icprobe1[-] ⟷ vscr1[+]
            t1[emitter] ⟷ reref[1]
            reref[2] ⟷ ieprobe1[+]
            veprobe1[-] ⟷ ieprobe1[-] ⟷ vscr1[-] ⟷ isrc1[-]
            t2 = bjt(typ, rb=rb, re=re, rc=rc)
            isrc2 = currentsource(ib)
            vscr2 = voltagesource(vce)
            veprobe2 = voltageprobe()
            vcprobe2 = voltageprobe()
            ieprobe2 = currentprobe()
            icprobe2 = currentprobe()
            t2[base] ⟷ isrc2[+] ⟷ veprobe2[+] ⟷ vcprobe2[+]
            t2[collector] ⟷ icprobe2[+]
            vcprobe2[-] ⟷ icprobe2[-] ⟷ vscr2[+]
            t2[emitter] ⟷ ieprobe2[+]
            veprobe2[-] ⟷ ieprobe2[-] ⟷ vscr2[-] ⟷ isrc2[-]
        end
        model = DiscreteModel(circ, 1)
        output = run!(model, zeros(0,1))
        @test output[1:4,:] ≈ output[5:8,:]
    end
end

@testset "MOSFET" begin
    for (typ, pol) in ((:n, 1), (:p, -1))
        circ = @circuit begin
            vgs = voltagesource(), [-] == gnd
            vds = voltagesource(), [-] == gnd
            J = mosfet(typ, vt=1, α=1e-4), [gate] == vgs[+], [drain] == vds[+]
            out = currentprobe(), [+] == J[source], [-] == gnd
        end
        model = DiscreteModel(circ, 1);
        y = run!(model, pol*[0 1 2 2 2; 5 5 0.5 1 1.5])
        @test y == pol*[0 0 1e-4*(1-0.5/2)*0.5 1e-4*(1-1/2)*1 1e-4/2*1^2]
    end
    for (typ, pol) in ((:n, 1), (:p, -1)), α in (1e-4, (0.0205, -0.0017)),
        vt in (1, (1.2078, 0.3238), (-1.2454, -0.199, -0.0483))
        circ = @circuit begin
            vgs = voltagesource(), [-] == gnd
            vds = voltagesource(), [-] == gnd
            J = mosfet(typ, vt=vt, α=α, λ=0.05), [gate] == vgs[+], [drain] == vds[+]
            out = currentprobe(), [+] == J[source], [-] == gnd
        end
        model = DiscreteModel(circ, 1);
        for vgs in range(0, stop=5, length=10), vds in range(0, stop=5, length=10)
            y = only(run!(model, pol*hcat([vgs; vds])))
            α´ = evalpoly(pol*vgs, (α...,))
            vt´ = evalpoly(pol*vgs, (vt...,))
            if vgs ≤ vt´
                @test y == 0
            elseif vds ≤ vgs - vt´
                @test y ≈ pol * α´ * (vgs - vt´ - vds / 2) * vds * (1 + 0.05 * vds)
            else
                @test y ≈ pol * α´ / 2 * (vgs - vt´)^2 * (1 + 0.05 * vds)
            end
        end
    end
end

@testset "op amp" begin
    for Amax in (10, Inf), GBP in (50e3, Inf)
        # test circuit: non-inverting amplifier in high shelving configuration
        circ = @circuit begin
            input = voltagesource(), [-] ⟷ gnd
            op = opamp(maxgain=Amax, gain_bw_prod=GBP), ["in+"] ⟷ input[+], ["out-"] ⟷ gnd
            r1 = resistor(109e3), [1] ⟷ op["out+"], [2] ⟷ op["in-"]
            r2 = resistor(1e3), [1] ⟷ op["in-"]
            c = capacitor(22e-9), [1] ⟷ r2[2], [2] ⟷ gnd
            output = voltageprobe(), [+] ⟷ op["out+"], [-] ⟷ gnd
        end
        model = DiscreteModel(circ, 1/44100)
        # obtain impulse response / transfer function
        u = [1; zeros(4095)]'
        y = run!(model, u)[1,:]
        Y = rfft(y)
        # inverse of op amp transfer function
        G⁻¹(s) = sqrt(1-1/Amax^2)*s/(2π*GBP) + 1/Amax
        #  feedback transfer function
        H(s) = (1e3*22e-9*s + 1) / ((109e3+1e3)*22e-9*s + 1)
        # overall transfer function evaluated taking frequency warping of
        # bilinear transform into account
        Yref = [let ω=2*44100*tan(π*k/length(y)); 1/(G⁻¹(im*ω) + H(im*ω)); end for k in eachindex(Y).-1]
        @test Y ≈ Yref
    end

    circ = @circuit begin
        input = voltagesource(), [-] ⟷ gnd
        op = opamp(Val{:macak}, 100, -3, 4), ["in+"] ⟷ input[+], ["in-"] ⟷ ["out-"] ⟷ gnd
        output = voltageprobe(), [+] ⟷ op["out+"], [-] ⟷ gnd
    end
    u = range(-1, stop=1, length=1000)
    model = DiscreteModel(circ, 1/44100)
    y = run!(model, u')[1,:]
    yref = [0.5*(4 + -3) + 0.5*(4 - -3)*tanh(100/(0.5*(4 - -3)) * u) for u in u]
    @test y ≈ yref
end

function checksteady!(model)
    x_steady = steadystate!(model)
    for s in model.solvers
        ACME.set_resabstol!(s, 1e-13)
    end
    run!(model, zeros(1, 1))
    return model.x ≈ x_steady
end

function linearization_error!(model, amplitude)
    linmodel = linearize(model)
    N = 50000
    u = [amplitude * sin(π/2 * n^2/N) for n in 0:N]'
    steadystate!(model)
    steadystate!(linmodel)
    y = run!(model, u)
    ylin = run!(linmodel, u)
    return maximum(abs, y-ylin)
end

@testset "examples" begin
    @testset "sallenkey" begin
        include(joinpath(dirname(@__FILE__), "..", "examples", "sallenkey.jl"))
         model=sallenkey()
         println("Running sallenkey")
         y = run!(model, sin.(2π*1000/44100*(0:44099)'); showprogress=false)
         @test size(y) == (1,44100)
         # TODO: further validate y
         @test checksteady!(model)
     end

     @testset "diodeclipper" begin
        include(joinpath(dirname(@__FILE__), "..", "examples", "diodeclipper.jl"))
        model=diodeclipper()
        println("Running diodeclipper")
        @test ACME.np(model, 1) == 1
        y = run!(model, sin.(2π*1000/44100*(0:44099)'); showprogress=false)
        @test size(y) == (1,44100)
        # TODO: further validate y
        @test checksteady!(model)

        @test linearization_error!(model, 1e-3) < 1e-15

        model = diodeclipper(solver=HomotopySolver{SimpleSolver})
        runner = ModelRunner(model, false)
        u = sin.(2π*1000/44100*(0:44099)')
        y = run!(runner, u)
        run!(runner, y, u)
        alloc = @allocated run!(runner, y, u)
        if alloc > 0
            warn("Allocated $alloc in run! of HomotopySolver{SimpleSolver} base ModelRunner")
        end
    end

    @testset "birdie" begin
        include(joinpath(dirname(@__FILE__), "..", "examples", "birdie.jl"))
        model=birdie(vol=0.8)
        ACME.solve(model.solvers[1], [0.003, -0.0002])
        @assert all(ACME.hasconverged, model.solvers)
        println("Running birdie with fixed vol")
        @test ACME.np(model, 1) == 2
        y = run!(model, sin.(2π*1000/44100*(0:44099)'); showprogress=false)
        @test size(y) == (1,44100)
        # TODO: further validate y
        @test checksteady!(model)

        @test linearization_error!(model, 1e-4) < 1e-7

        model=birdie()
        println("Running birdie with varying vol")
        @test ACME.np(model, 1) == 3
        y = run!(model, [sin.(2π*1000/44100*(0:44099)'); range(1, stop=0, length=44100)']; showprogress=false)
        @test size(y) == (1,44100)
        # TODO: further validate y
    end

    @testset "superover" begin
        include(joinpath(dirname(@__FILE__), "..", "examples", "superover.jl"))
        model=superover(drive=1.0, tone=1.0, level=1.0)
        println("Running superover with fixed potentiometer values")
        @test ACME.np(model, 1) == 5
        y = run!(model, sin.(2π*1000/44100*(0:44099)'); showprogress=false)
        @test size(y) == (1,44100)
        # TODO: further validate y
        @test checksteady!(model)
        @test linearization_error!(model, 1e-4) < 1e-4 # SuperOver really is not very linear...

        circ=superover(Circuit, drive=1.0, tone=1.0, level=1.0)
        println("Running simplified superover with fixed potentiometer values")
        add!(circ, :vbsrc, voltagesource(4.5))
        connect!(circ, (:vbsrc, :+), :vb)
        connect!(circ, (:vbsrc, :-), :gnd)
        model = DiscreteModel(circ, 1/44100)
        @test ACME.np(model, 1) == 2
        @test ACME.np(model, 2) == 1
        @test ACME.np(model, 3) == 2
        y = run!(model, sin.(2π*1000/44100*(0:44099)'); showprogress=false)
        @test size(y) == (1,44100)
        # TODO: further validate y
        @test_broken checksteady!(model)
        @test_broken linearization_error!(model, 1e-4) < 1e-4 # SuperOver really is not very linear...

        println("Running simplified, non-decomposed superover with fixed potentiometer values")
        model = DiscreteModel(circ, 1/44100, decompose_nonlinearity=false)
        @test ACME.np(model, 1) == 5
        y = run!(model, sin.(2π*1000/44100*(0:44099)'); showprogress=false)
        @test size(y) == (1,44100)
        # TODO: further validate y
        @test checksteady!(model)
        @test linearization_error!(model, 1e-4) < 1e-4 # SuperOver really is not very linear...

        model=superover()
        println("Running superover with varying potentiometer values")
        @test ACME.np(model, 1) == 11
        y = run!(model, [sin.(2π*1000/44100*(0:999)'); range(1, stop=0, length=1000)'; range(0, stop=1, length=1000)'; range(1, stop=0, length=1000)']; showprogress=false)
        @test size(y) == (1,1000)
        # TODO: further validate y

        circ=superover(Circuit)
        println("Running simplified superover with varying potentiometer values")
        add!(circ, :vbsrc, voltagesource(4.5))
        connect!(circ, (:vbsrc, :+), :vb)
        connect!(circ, (:vbsrc, :-), :gnd)
        model = DiscreteModel(circ, 1/44100)
        @test ACME.np(model, 1) == 2
        @test ACME.np(model, 2) == 2
        @test ACME.np(model, 3) == 2
        @test ACME.np(model, 4) == 4
        y = run!(model, [sin.(2π*1000/44100*(0:999)'); range(1, stop=0, length=1000)'; range(0, stop=1, length=1000)'; range(1, stop=0, length=1000)']; showprogress=false)
        @test size(y) == (1,1000)
        # TODO: further validate y
    end
end
