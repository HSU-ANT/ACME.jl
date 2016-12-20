# Copyright 2015, 2016 Martin Holters
# See accompanying license file.

using ACME
using Base.Test

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
    @test_approx_eq A*y x
    copy!(y, x)
    ACME.solve!(solver, y, y)
    @test_approx_eq A*y x
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
    @test !isempty(search(convert(ASCIIString, readavailable(rd)), "WARNING"))
    redirect_stderr(orig_stderr)
end

let circ = Circuit(), d = diode(), src=currentsource()
    connect!(circ, d[:+], src[:+])
    connect!(circ, d[:-], src[:-])
    model = DiscreteModel(circ, 1)
    @test_throws ErrorException run!(model, hcat([Inf]))
    orig_stderr = STDERR
    rd, wr = redirect_stderr()
    @test run!(model, hcat([-1.0])) == zeros(0, 1)
    # should warn because solution exists as diode cannot reach reverse current of 1A
    @test !isempty(search(convert(ASCIIString, readavailable(rd)), "WARNING"))
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
    @test_approx_eq sum(abs2, p - best_p) sum(abs2, p - ps[:, idx])
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
    @test_approx_eq_eps y[1] v_d 1e-6
end

function checksteady!(model)
    x_steady = steadystate!(model)
    ACME.set_resabs2tol!(model.solver, 1e-25)
    run!(model, zeros(1, 1))
    @test_approx_eq model.x x_steady
end

include("../examples/sallenkey.jl")
let model=sallenkey()
    println("Running sallenkey")
    y = run!(model, map(sin, 2π*1000/44100*(0:44099)'))
    @test size(y) == (1,44100)
    # TODO: further validate y

    # cannot check steady state: steadystate() does not work for matrix A having
    # eigenvalue 1
end

include("../examples/diodeclipper.jl")
let model=diodeclipper()
    println("Running diodeclipper")
    @test ACME.np(model) == 1
    y = run!(model, map(sin, 2π*1000/44100*(0:44099)'))
    @test size(y) == (1,44100)
    # TODO: further validate y
    checksteady!(model)
end

include("../examples/birdie.jl")
let model=birdie(vol=0.8)
    ACME.solve(model.solver, [0.003, -0.0002])
    @assert ACME.hasconverged(model.solver)
    println("Running birdie with fixed vol")
    @test ACME.np(model) == 2
    y = run!(model, map(sin, 2π*1000/44100*(0:44099)'))
    @test size(y) == (1,44100)
    # TODO: further validate y
    checksteady!(model)
end
let model=birdie()
    println("Running birdie with varying vol")
    @test ACME.np(model) == 3
    y = run!(model, [map(sin, 2π*1000/44100*(0:44099).'); linspace(1,0,44100).'])
    @test size(y) == (1,44100)
    # TODO: further validate y
end

include("../examples/superover.jl")
let model=superover(drive=1.0, tone=1.0, level=1.0)
    println("Running superover with fixed potentiometer values")
    @test ACME.np(model) == 5
    y = run!(model, map(sin, 2π*1000/44100*(0:44099)'))
    @test size(y) == (1,44100)
    # TODO: further validate y
    checksteady!(model)
end
let model=superover()
    println("Running superover with varying potentiometer values")
    @test ACME.np(model) == 11
    y = run!(model, [map(sin, 2π*1000/44100*(0:999)'); linspace(1,0,1000).'; linspace(0,1,1000).'; linspace(1,0,1000).'])
    @test size(y) == (1,1000)
    # TODO: further validate y
end
