module ACME

type Element
  mv :: SparseMatrixCSC{Number,Int}
  mi :: SparseMatrixCSC{Number,Int}
  mx :: SparseMatrixCSC{Number,Int}
  mxd :: SparseMatrixCSC{Number,Int}
  mq :: SparseMatrixCSC{Number,Int}
  mu :: SparseMatrixCSC{Number,Int}
  u0 :: SparseMatrixCSC{Number,Int}
  pv :: SparseMatrixCSC{Number,Int}
  pi :: SparseMatrixCSC{Number,Int}
  px :: SparseMatrixCSC{Number,Int}
  pxd :: SparseMatrixCSC{Number,Int}
  pq :: SparseMatrixCSC{Number,Int}
  nonlinear_eq :: Expr
  pins :: Dict{Symbol, Vector{(Int, Int)}}

  function Element(;args...)
    sizes = (Symbol=>Int)[:n0 => 1]

    function update_sizes(mat, syms)
      for (sym, s) in zip(syms, size(mat))
        if !haskey(sizes, sym)
          sizes[sym] = s
        elseif sizes[sym] ≠ s
          error("Inconsistent sizes for ", sym)
        end
      end
    end

    function make_pin_dict(syms)
      dict = (Symbol=>Vector{(Int, Int)})[]
      for i in 1:length(syms)
        branch = div(i+1, 2)
        polarity = 2mod(i, 2) - 1
        push!(get!(dict, symbol(syms[i]), []), (branch, polarity))
      end
      dict
    end
    make_pin_dict(dict::Dict) = dict

    const mat_dims = [ :mv => (:nl,:nb), :mi => (:nl,:nb), :mx => (:nl,:nx),
                       :mxd => (:nl,:nx), :mq => (:nl,:nq), :mu => (:nl,:nu),
                       :u0 => (:nl, :n0),
                       :pv => (:ny,:nb), :pi => (:ny,:nb), :px => (:ny,:nx),
                       :pxd => (:ny,:nx), :pq => (:ny,:nq) ]

    elem = new()
    for (key, val) in args
      if haskey(mat_dims, key)
        val = sparse([val])
        update_sizes (val, mat_dims[key])
      elseif key == :pins
        val = make_pin_dict (val)
      end
      elem.(key) = val
    end
    for (m, ns) in mat_dims
      if !isdefined(elem, m)
        elem.(m) = spzeros(Int, get(sizes, ns[1], 0), get(sizes, ns[2], 0))
      end
    end
    if !isdefined(elem, :nonlinear_eq)
      elem.nonlinear_eq = Expr(:block)
    end
    if !isdefined(elem, :pins)
      elem.pins = make_pin_dict(map(string,1:2nb(elem)))
    end
    elem
  end
end

for (n,m) in [:nb => :mv, :nx => :mx, :nq => :mq, :nu => :mu]
  @eval ($n)(e::Element) = size(e.$m)[2]
end
nl(e::Element) = size(e.mv)[1]
ny(e::Element) = size(e.pv)[1]
nn(e::Element) = nb(e) + nx(e) + nq(e) - nl(e)

include("elements.jl")

end # module
