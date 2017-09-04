# Copyright 2015, 2016, 2017 Martin Holters
# See accompanying license file.

using ACME
using Base.Test
using Compat
using ProgressMeter

tv, ti = ACME.topomat(sparse([1 -1 1; -1 1 -1]))
@test tv*ti'==spzeros(2,1)

# Pathological cases for topomat:
# two nodes, one loop branch (short-circuited) -> voltage==0, current arbitrary
@test ACME.topomat(spzeros(Int, 2, 1)) == (speye(1), spzeros(0, 1))
# two nodes, one branch between them -> voltage arbitrary, current==0
@test ACME.topomat(sparse([1,2], [1,1], [1,-1])) == (spzeros(0, 1), speye(1))

let solver = ACME.LinearSolver(3)
    A = [1.0 0.5 0.4; 2.0 4.0 1.7; 4.0 7.0 9.1]
    @test ACME.setlhs!(solver, A)
    x = rand(3)
    y = similar(x)
    ACME.solve!(solver, y, x)
    @test A*y ≈ x
    copy!(y, x)
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

let circ = Circuit()
    model=DiscreteModel(circ, 1)
    @test run!(model, zeros(0, 20)) == zeros(0, 20)
end

let circ = Circuit(), r = resistor(0)
    connect!(circ, r[1], r[2])
    model = DiscreteModel(circ, 1)
    @test run!(model, zeros(0, 20)) == zeros(0, 20)
end

let circ = Circuit(), r = resistor(0), probe = currentprobe()
    connect!(circ, r[1], probe[:+])
    connect!(circ, r[2], probe[:-])
    orig_stderr = STDERR
    rd, wr = redirect_stderr()
    model = DiscreteModel(circ, 1)
    # should warn because output is indeterminate
    @test !isempty(search(convert(Compat.ASCIIString, readavailable(rd)), "WARNING"))
    redirect_stderr(orig_stderr)
end

let circ = Circuit(), d = diode(), src=currentsource(), probe=voltageprobe()
    connect!(circ, d[:+], src[:+])
    connect!(circ, d[:-], src[:-])
    connect!(circ, d[:+], probe[:+])
    connect!(circ, d[:-], probe[:-])
    model = DiscreteModel(circ, 1)
    @test ACME.nn(model) == 1
    y = run!(model, [1.0 1.0])
    @test size(y) == (1, 2)
    @test y[1,1] == y[1,2]
    @test_throws ErrorException run!(model, hcat([Inf]))
    orig_stderr = STDERR
    rd, wr = redirect_stderr()
    @test size(run!(model, hcat([-1.0]))) == (1, 1)
    # should warn because solution exists as diode cannot reach reverse current of 1A
    @test !isempty(search(convert(Compat.ASCIIString, readavailable(rd)), "WARNING"))
    redirect_stderr(orig_stderr)
end

for num = 1:50
    let ps = rand(4, num)
        t = ACME.KDTree(ps)
        for i in 1:size(ps)[2]
            idx = ACME.indnearest(t, ps[:,i])
            @test ps[:,i] == ps[:,idx]
        end
    end
end

let ps = rand(6, 10000)
    t = ACME.KDTree(ps)
    p = rand(6)
    best_p = ps[:,indmin(sum(abs2, broadcast(-, ps, p),1))]
    idx = ACME.indnearest(t, p)
    @test sum(abs2, p - best_p) ≈ sum(abs2, p - ps[:, idx])
end

let nleq = ACME.ParametricNonLinEq((res, J, scratch, z) ->
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

let a = Rational{BigInt}[1 1 1; 1 1 2; 1 2 1; 1 2 2; 2 1 1; 2 1 2],
    b = Rational{BigInt}[1 2 3 4 5 6; 6 5 4 3 2 1; 1 0 1 0 1 0]
    nullspace = ACME.gensolve(sparse(a'), spzeros(size(a, 2), 0))[2]
    @test nullspace'*a == spzeros(3, 3)
    c, f = ACME.rank_factorize(sparse(a * b))
    @test c*f == a*b
end

let a = Rational{BigInt}[1 1 1; 1 1 2; 1 2 1; 1 2 2; 2 1 1; 2 1 2],
    b = Rational{BigInt}[1 2 3 4 5 6; 6 5 4 3 2 1; 1 0 1 0 1 0],
    fq = Rational{BigInt}[1 0 0; 10 0 0; 0 1 0; 0 10 0; 0 0 1; 0 0 10],
    z = Rational{BigInt}[1 2 0 0 2 1; 0 1 2 2 0 1; 0 0 1 0 1 1]
    mats = Dict{Symbol,Array}(:dq_full => a * b, :eq_full => zeros(Rational{BigInt},6,0), :fq => fq)
    mats[:dq_fulls]=Matrix[mats[:dq_full]]
    mats[:eq_fulls]=Matrix[mats[:eq_full]]
    mats[:fqprev_fulls]=Matrix[mats[:eq_full]]
    mats[:fqs]=Matrix[mats[:fq]]
    ACME.reduce_pdims!(mats)
    @test size(mats[:pexps][1], 2) == 3
    @test mats[:pexps][1] * mats[:dqs][1] == mats[:dq_fulls][1]
    mats = Dict{Symbol,Array}(:dq_full => a * b + fq * z, :eq_full => zeros(Rational{BigInt},6,0), :fq => fq)
    mats[:dq_fulls]=Matrix[mats[:dq_full]]
    mats[:eq_fulls]=Matrix[mats[:eq_full]]
    mats[:fqprev_fulls]=Matrix[mats[:eq_full]]
    mats[:fqs]=Matrix[mats[:fq]]
    ACME.reduce_pdims!(mats)
    @test size(mats[:pexps][1], 2) == 3
end

let circ = Circuit()
    src1 = voltagesource()
    probe1 = currentprobe()
    d1 = diode()
    d2 = diode()
    connect!(circ, src1[:+], d1[:+])
    connect!(circ, d1[:-], d2[:+])
    connect!(circ, d2[:-], probe1[:+])
    connect!(circ, probe1[:-], src1[:-])
    src2 = voltagesource()
    probe2 = currentprobe()
    d3 = diode()
    connect!(circ, src2[:+], d3[:+])
    connect!(circ, d3[:-], probe2[:+])
    connect!(circ, probe2[:-], src2[:-])
    model = DiscreteModel(circ, 1)
    y = run!(model, hcat([2.0; 1.0]))
    # single diode is extracted first, although it was added last
    @test ACME.nn(model, 1) == 1
    @test ACME.nn(model, 2) == 2
    @test y[1] ≈ y[2] ≈ 1e-12*(exp(1/25e-3)-1)
end

# sources and probes with internal resistance/conductance
let circ = Circuit(), src=currentsource(100e-3, gp=1//100000), probe=voltageprobe()
    connect!(circ, src[:+], probe[:+])
    connect!(circ, src[:-], probe[:-])
    model = DiscreteModel(circ, 1)
    @test run!(model, zeros(0,1)) ≈ [100000*100e-3]
end
let circ = Circuit(), src=currentsource(gp=1//100000), probe=voltageprobe()
    connect!(circ, src[:+], probe[:+])
    connect!(circ, src[:-], probe[:-])
    model = DiscreteModel(circ, 1)
    @test run!(model, hcat([100e-3])) ≈ [100000*100e-3]
end
let circ = Circuit(), src=currentsource(100e-3), probe=voltageprobe(gp=1//100000)
    connect!(circ, src[:+], probe[:+])
    connect!(circ, src[:-], probe[:-])
    model = DiscreteModel(circ, 1)
    @test run!(model, zeros(0,1)) ≈ [100000*100e-3]
end
let circ = Circuit(), src=voltagesource(10, rs=100000), probe=currentprobe()
    connect!(circ, src[:+], probe[:+])
    connect!(circ, src[:-], probe[:-])
    model = DiscreteModel(circ, 1)
    @test run!(model, zeros(0,1)) ≈ [10/100000]
end
let circ = Circuit(), src=voltagesource(rs=100000), probe=currentprobe()
    connect!(circ, src[:+], probe[:+])
    connect!(circ, src[:-], probe[:-])
    model = DiscreteModel(circ, 1)
    @test run!(model, hcat([10.0])) ≈ [10/100000]
end
let circ = Circuit(), src=voltagesource(10), probe=currentprobe(rs=100000)
    connect!(circ, src[:+], probe[:+])
    connect!(circ, src[:-], probe[:-])
    model = DiscreteModel(circ, 1)
    @test run!(model, zeros(0,1)) ≈ [10/100000]
end

# BJT Ebers-Moll model
let isc=1e-6, ise=2e-6, ηc=1.1, ηe=1.0, βf=100, βr=10
    for (typ, ib) in ((:npn, 1e-3), (:pnp, -1e-3))
        t = bjt(typ, isc=isc, ise=ise, ηc=ηc, ηe=ηe, βf=βf, βr=βr)
        isrc = currentsource()
        vsrc = voltagesource()
        veprobe = voltageprobe()
        vcprobe = voltageprobe()
        ieprobe = currentprobe()
        icprobe = currentprobe()
        circ = Circuit()
        add!(circ, veprobe, vcprobe, ieprobe, icprobe)
        connect!(circ, t[:base], isrc[:+], veprobe[:+], vcprobe[:+])
        connect!(circ, t[:collector], icprobe[:+])
        connect!(circ, vcprobe[:-], icprobe[:-], vsrc[:+])
        connect!(circ, t[:emitter], ieprobe[:+])
        connect!(circ, veprobe[:-], ieprobe[:-], vsrc[:-], isrc[:-])
        model = DiscreteModel(circ, 1)
        N = 100
        output = run!(model, [linspace(0, ib, N).'; linspace(1, -1, N÷2).' linspace(-1, 1, N÷2).'])
        if typ == :pnp
            output = -output
        end
        for n in 1:N
            ve, vc, ie, ic = output[:,n]
            @test isapprox(ie, ise*(exp(ve/(ηe*25e-3))-1) - βr/(1+βr)*isc*(exp(vc/(ηc*25e-3))-1), atol=1e-10)
            @test isapprox(ic, -βf/(1+βf)*ise*(exp(ve/(ηe*25e-3))-1) + isc*(exp(vc/(ηc*25e-3))-1), atol=1e-10)
        end
    end
end
# BJT Gummel-Poon model
let isc=1e-6, ise=2e-6, ηc=1.1, ηe=1.0, βf=100, βr=10, ηcl=1.2, ηel=1.3
    prog = Progress(2^9, "Testing Gummel-Poon model: ")
    for ile in (0, 50e-9), ilc in (0, 100e-9),
            ηcl in (ηc, 1.2), ηel in (ηe, 1.1),
            vaf in (Inf, 10), var in (Inf, 50),
            ikf in (Inf, 50e-3), ikr in (Inf, 500e-3),
            (typ, ib) in ((:npn, 1e-3), (:pnp, -1e-3))
        t = bjt(typ, isc=isc, ise=ise, ηc=ηc, ηe=ηe, βf=βf, βr=βr, ile=ile,
                ilc=ilc, ηcl=ηcl, ηel=ηel, vaf=vaf, var=var, ikf=ikf, ikr=ikr)
        isrc = currentsource()
        vsrc = voltagesource()
        veprobe = voltageprobe()
        vcprobe = voltageprobe()
        ieprobe = currentprobe()
        icprobe = currentprobe()
        circ = Circuit()
        add!(circ, veprobe, vcprobe, ieprobe, icprobe)
        connect!(circ, t[:base], isrc[:+], veprobe[:+], vcprobe[:+])
        connect!(circ, t[:collector], icprobe[:+])
        connect!(circ, vcprobe[:-], icprobe[:-], vsrc[:+])
        connect!(circ, t[:emitter], ieprobe[:+])
        connect!(circ, veprobe[:-], ieprobe[:-], vsrc[:-], isrc[:-])
        model = DiscreteModel(circ, 1)
        N = 100
        output = run!(model, [linspace(0, ib, N).'; linspace(1, -1, N÷2).' linspace(-1, 1, N÷2).'])
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
end
# compare internal to external terminal resistances
let rb=100, re=10, rc=20
    for (typ, ib, vce) in ((:npn, 1e-3, 1), (:pnp, -1e-3, -1))
        t1 = bjt(typ)
        rbref=resistor(rb)
        rcref=resistor(rc)
        reref=resistor(re)
        isrc1 = currentsource(ib)
        vscr1 = voltagesource(vce)
        veprobe1 = voltageprobe()
        vcprobe1 = voltageprobe()
        ieprobe1 = currentprobe()
        icprobe1 = currentprobe()
        t2 = bjt(typ, rb=rb, re=re, rc=rc)
        isrc2 = currentsource(ib)
        vscr2 = voltagesource(vce)
        veprobe2 = voltageprobe()
        vcprobe2 = voltageprobe()
        ieprobe2 = currentprobe()
        icprobe2 = currentprobe()
        circ = Circuit()
        add!(circ, veprobe1, vcprobe1, ieprobe1, icprobe1)
        connect!(circ, t1[:base], rbref[1])
        connect!(circ, rbref[2], isrc1[:+], veprobe1[:+], vcprobe1[:+])
        connect!(circ, t1[:collector], rcref[1])
        connect!(circ, rcref[2], icprobe1[:+])
        connect!(circ, vcprobe1[:-], icprobe1[:-], vscr1[:+])
        connect!(circ, t1[:emitter], reref[1])
        connect!(circ, reref[2], ieprobe1[:+])
        connect!(circ, veprobe1[:-], ieprobe1[:-], vscr1[:-], isrc1[:-])
        add!(circ, veprobe2, vcprobe2, ieprobe2, icprobe2)
        connect!(circ, t2[:base], isrc2[:+], veprobe2[:+], vcprobe2[:+])
        connect!(circ, t2[:collector], icprobe2[:+])
        connect!(circ, vcprobe2[:-], icprobe2[:-], vscr2[:+])
        connect!(circ, t2[:emitter], ieprobe2[:+])
        connect!(circ, veprobe2[:-], ieprobe2[:-], vscr2[:-], isrc2[:-])
        model = DiscreteModel(circ, 1)
        output = run!(model, zeros(0,1))
        @test output[1:4,:] ≈ output[5:8,:]
    end
end

# simple circuit: resistor and diode in series, driven by constant voltage,
# chosen such that a prescribe current flows
let i = 1e-3, r=10e3, is=1e-12
    v_r = i*r
    v_d = 25e-3 * log(i/is+1)
    vsrc = voltagesource(v_r + v_d)
    r1 = resistor(r)
    d = diode(is=is)
    vprobe = voltageprobe()
    circ = Circuit()
    add!(circ, vsrc)
    add!(circ, r1, d)
    connect!(circ, vsrc[:+], :vcc)
    connect!(circ, vsrc[:-], :gnd)
    connect!(circ, r1[1], :vcc)
    connect!(circ, d[:-], vprobe[:-], :gnd)
    connect!(circ, r1[2], d[:+], vprobe[:+])
    model = DiscreteModel(circ, 1)
    y = run!(model, zeros(0, 1))
    @test y[1] ≈ v_d
end

function checksteady!(model)
    x_steady = steadystate!(model)
    for s in model.solvers
        ACME.set_resabstol!(s, 1e-13)
    end
    run!(model, zeros(1, 1))
    @test model.x ≈ x_steady
end

include("../examples/sallenkey.jl")
let model=sallenkey()
    println("Running sallenkey")
    y = run!(model, map(sin, 2π*1000/44100*(0:44099)'); showprogress=false)
    @test size(y) == (1,44100)
    # TODO: further validate y

    # cannot check steady state: steadystate() does not work for matrix A having
    # eigenvalue 1
end

include("../examples/diodeclipper.jl")
let model=diodeclipper()
    println("Running diodeclipper")
    @test ACME.np(model, 1) == 1
    y = run!(model, map(sin, 2π*1000/44100*(0:44099)'); showprogress=false)
    @test size(y) == (1,44100)
    # TODO: further validate y
    checksteady!(model)

    linmodel = linearize(model)
    N = 50000
    u = [1e-3 * sin(π/2 * n^2/N) for n in 0:N]'
    steadystate!(model)
    steadystate!(linmodel)
    y = run!(model, u)
    ylin = run!(linmodel, u)
    @test y ≈ ylin
end
let circ = diodeclipper(Circuit)
    model = DiscreteModel(circ, 44100, ACME.HomotopySolver{ACME.SimpleSolver})
    runner = ModelRunner(model, false)
    u = map(sin, 2π*1000/44100*(0:44099)')
    y = run!(runner, u)
    run!(runner, y, u)
    alloc = @allocated run!(runner, y, u)
    if alloc > 0
        warn("Allocated $alloc in run! of HomotopySolver{SimpleSolver} base ModelRunner")
    end
end

include("../examples/birdie.jl")
let model=birdie(vol=0.8)
    ACME.solve(model.solvers[1], [0.003, -0.0002])
    @assert all(ACME.hasconverged, model.solvers)
    println("Running birdie with fixed vol")
    @test ACME.np(model, 1) == 2
    y = run!(model, map(sin, 2π*1000/44100*(0:44099)'); showprogress=false)
    @test size(y) == (1,44100)
    # TODO: further validate y
    checksteady!(model)

    linmodel = linearize(model)
    N = 50000
    u = [1e-4 * sin(π/2 * n^2/N) for n in 0:N]'
    steadystate!(model)
    steadystate!(linmodel)
    y = run!(model, u)
    ylin = run!(linmodel, u)
    @test maximum(abs, y-ylin) < 1e-6
end
let model=birdie()
    println("Running birdie with varying vol")
    @test ACME.np(model, 1) == 3
    y = run!(model, [map(sin, 2π*1000/44100*(0:44099).'); linspace(1,0,44100).']; showprogress=false)
    @test size(y) == (1,44100)
    # TODO: further validate y
end

include("../examples/dasplus.jl")
let model=dasplus(gain=0.8, level=0.5)
    println("Running dasplus with fixed potentiometer values")
    @test ACME.np(model, 1) == 1
    @test ACME.np(model, 2) == 1
    y = run!(model, map(sin, 2π*1000/44100*(0:44099)'); showprogress=false)
    @test size(y) == (1,44100)
    # TODO: further validate y
    checksteady!(model)
end
let model=dasplus()
    println("Running dasplus with varying potentiometer values")
    @test ACME.np(model, 1) == 3
    @test ACME.np(model, 2) == 2
    y = run!(model, [map(sin, 2π*1000/44100*(0:44099).'); linspace(1,0,44100).'; linspace(0,1,44100).']; showprogress=false)
    @test size(y) == (1,44100)
    # TODO: further validate y
end

include("../examples/superover.jl")
let model=superover(drive=1.0, tone=1.0, level=1.0)
    println("Running superover with fixed potentiometer values")
    @test ACME.np(model, 1) == 5
    y = run!(model, map(sin, 2π*1000/44100*(0:44099)'); showprogress=false)
    @test size(y) == (1,44100)
    # TODO: further validate y
    checksteady!(model)

    linmodel = linearize(model)
    N = 50000
    u = [1e-4 * sin(π/2 * n^2/N) for n in 0:N]'
    steadystate!(model)
    steadystate!(linmodel)
    y = run!(model, u)
    ylin = run!(linmodel, u)
    @test maximum(abs, y-ylin) < 2e-4 # SuperOver really is not very linear...
end
let circ=superover(Circuit, drive=1.0, tone=1.0, level=1.0)
    println("Running simplified superover with fixed potentiometer values")
    vbsrc = voltagesource(4.5)
    connect!(circ, vbsrc[:+], :vb)
    connect!(circ, vbsrc[:-], :gnd)
    model = DiscreteModel(circ, 1/44100)
    @test ACME.np(model, 1) == 2
    @test ACME.np(model, 2) == 1
    @test ACME.np(model, 3) == 2
    y = run!(model, map(sin, 2π*1000/44100*(0:44099)'); showprogress=false)
    @test size(y) == (1,44100)
    # TODO: further validate y
    checksteady!(model)

    linmodel = linearize(model)
    N = 50000
    u = [1e-4 * sin(π/2 * n^2/N) for n in 0:N]'
    steadystate!(model)
    steadystate!(linmodel)
    y = run!(model, u)
    ylin = run!(linmodel, u)
    @test maximum(abs, y-ylin) < 2e-4 # SuperOver really is not very linear...

    println("Running simplified, non-decomposed superover with fixed potentiometer values")
    model = DiscreteModel(circ, 1/44100, decompose_nonlinearity=false)
    @test ACME.np(model, 1) == 5
    y = run!(model, map(sin, 2π*1000/44100*(0:44099)'); showprogress=false)
    @test size(y) == (1,44100)
    # TODO: further validate y
    checksteady!(model)

    linmodel = linearize(model)
    N = 50000
    u = [1e-4 * sin(π/2 * n^2/N) for n in 0:N]'
    steadystate!(model)
    steadystate!(linmodel)
    y = run!(model, u)
    ylin = run!(linmodel, u)
    @test maximum(abs, y-ylin) < 2e-4 # SuperOver really is not very linear...
end
let model=superover()
    println("Running superover with varying potentiometer values")
    @test ACME.np(model, 1) == 11
    y = run!(model, [map(sin, 2π*1000/44100*(0:999)'); linspace(1,0,1000).'; linspace(0,1,1000).'; linspace(1,0,1000).']; showprogress=false)
    @test size(y) == (1,1000)
    # TODO: further validate y
end
let circ=superover(Circuit)
    println("Running simplified superover with varying potentiometer values")
    vbsrc = voltagesource(4.5)
    connect!(circ, vbsrc[:+], :vb)
    connect!(circ, vbsrc[:-], :gnd)
    model = DiscreteModel(circ, 1/44100)
    @test ACME.np(model, 1) == 2
    @test ACME.np(model, 2) == 2
    @test ACME.np(model, 3) == 2
    @test ACME.np(model, 4) == 4
    y = run!(model, [map(sin, 2π*1000/44100*(0:999)'); linspace(1,0,1000).'; linspace(0,1,1000).'; linspace(1,0,1000).']; showprogress=false)
    @test size(y) == (1,1000)
    # TODO: further validate y
end
