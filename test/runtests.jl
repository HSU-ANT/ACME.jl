# Copyright 2015, 2016 Martin Holters
# See accompanying license file.

using ACME
using Base.Test

tv, ti = ACME.topomat(sparse([1 -1 1; -1 1 -1]))
@test tv*ti'==spzeros(2,1)

# Pathological cases for topomat:
# two nodes, one loop branch (short-circuited) -> voltage==0, current arbitrary
@test ACME.topomat(spzeros(Int, 2, 1)) == (sparse([1]), sparse([]))
# two nodes, one branch between them -> voltage arbitrary, current==0
@test ACME.topomat(sparse([1,2], [1,1], [1,-1])) == (sparse([]), sparse([1]))

let circ = Circuit()
    model=DiscreteModel(circ, 1.)
    @test run(model, zeros(0, 20)) == zeros(0, 20)
end

for num = 1:50
    let ps = rand(4, num)
        t = ACME.KDTree(ps)
        for i in 1:size(ps)[2]
            idx = ACME.indnearest(t, ps[:,i])[1]
            @test ps[:,i] == ps[:,idx]
        end
    end
end

let ps = rand(6, 10000)
    t = ACME.KDTree(ps)
    p = rand(6)
    best_p = ps[:,indmin(sumabs2(broadcast(-, ps, p),1))]
    idx = ACME.indnearest(t, p)[1]
    @test_approx_eq sumabs2(p - best_p) sumabs2(p - ps[:, idx])
end

# simple circuit: resistor and diode in series, driven by constant voltage,
# chosen such that a prescribe current flows
let i = 1e-3, r=10e3, is=1e-12
    v_r = i*r
    v_d = 25e-3 * log(i/is+1)
    vsrc = voltagesource(v_r + v_d)
    r1 = resistor(r)
    d = diode(is)
    vprobe = voltageprobe()
    circ = Circuit()
    connect!(circ, vsrc[:+], :vcc)
    connect!(circ, vsrc[:-], :gnd)
    connect!(circ, r1[1], :vcc)
    connect!(circ, d[:-], vprobe[:-], :gnd)
    connect!(circ, r1[2], d[:+], vprobe[:+])
    model = DiscreteModel(circ, 1.)
    y = run(model, zeros(0, 1))
    @test_approx_eq_eps y[1] v_d 1e-6
end
