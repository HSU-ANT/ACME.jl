# Copyright 2015, 2016, 2017, 2018 Martin Holters
# See accompanying license file.

include("checklic.jl")

using ACME
using Compat
using Compat.Test
using ProgressMeter

if VERSION ≥ v"0.7.0-DEV.3389"
    using SparseArrays
end
if VERSION < v"0.7.0-DEV.3986"
    range(start; stop=error("missing stop"), length=error("missing length")) = linspace(start, stop, length)
end
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
            vsrc = voltagesource(v_r + v_d), [+] ⟷ vcc,[-] ⟷ gnd
            r1 = resistor(r)
            d = diode(is=is), [-] ⟷ gnd, [+] ⟷ r1[2]
            vprobe = voltageprobe(), [-] ⟷ gnd, [+] ⟷ r1[2]
            r1[1] ⟷ vcc
        end
        model = DiscreteModel(circ, 1)
        y = run!(model, zeros(0, 1))
        @test y[1] ≈ v_d
    end
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
        @static if VERSION ≥ v"0.7.0-DEV.2988"
            @test_logs (:warn, "Model output depends on indeterminate quantity") DiscreteModel(circ, 1)
        else
            @static if VERSION ≥ v"0.6.0"
                @test_warn "output depends on indeterminate quantity" DiscreteModel(circ, 1)
            else
                orig_stderr = STDERR
                rd, wr = redirect_stderr()
                DiscreteModel(circ, 1)
                # should warn because output is indeterminate
                @test !isempty(search(String(readavailable(rd)), "WARNING"))
                redirect_stderr(orig_stderr)
            end
        end
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
        @static if VERSION ≥ v"0.7.0-DEV.2988"
            @test(size(@test_logs((:warn, "Failed to converge while solving non-linear equation."), run!(model, hcat([-1.0])))) == (1, 1))
        else
            @static if VERSION ≥ v"0.6.0"
                @test_warn("Failed to converge", @test size(run!(model, hcat([-1.0]))) == (1, 1))
            else
                orig_stderr = STDERR
                rd, wr = redirect_stderr()
                @test size(run!(model, hcat([-1.0]))) == (1, 1)
                # should warn because solution exists as diode cannot reach reverse current of 1A
                @test !isempty(search(String(readavailable(rd)), "WARNING"))
                redirect_stderr(orig_stderr)
            end
        end
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
        if isdefined(@__MODULE__, :argmin) # since 0.7.0-DEV.3516
            best_p = ps[:,argmin(vec(sum(abs2, ps .- p, 1)))]
        else
            best_p = ps[:,indmin(vec(sum(abs2, ps .- p, 1)))]
        end
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
    nullspace = ACME.gensolve(sparse(a'), spzeros(size(a, 2), 0))[2]
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
