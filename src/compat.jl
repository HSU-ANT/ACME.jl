using Compat

if VERSION ≥ v"0.6.0"
    macro pfunction(sig, params, body)
        esc(Expr(:function, Expr(:where, sig, params.args...), body))
    end
else
    macro pfunction(sig, params, body)
        esc(Expr(:function,
                 Expr(:call, Expr(:curly, sig.args[1], params.args...),
                      sig.args[2:end]...),
                 body))
    end
end

if VERSION ≥ v"0.6.0"
    @eval macro $(:struct)(head, body)
        Expr(Meta.parse("struct Foo end").head, false, esc(head), Expr(:block, esc.(body.args)...))
    end
    macro mutable_struct(head, body)
        Expr(Meta.parse("struct Foo end").head, true, esc(head), Expr(:block, esc.(body.args)...))
    end
else
    @eval macro $(:struct)(head, body)
        Expr(:type, false, esc(head), Expr(:block, esc.(body.args)...))
    end
    macro mutable_struct(head, body)
        Expr(:type, true, esc(head), Expr(:block, esc.(body.args)...))
    end
end

if !isdefined(@__MODULE__, :copyto!) # prior to 0.7.0-DEV.3057
    global const copyto! = Base.copy!
end

if !isdefined(@__MODULE__, :findall) # prior to 0.7.0-DEV.3415
    global const findall = Base.find
end

if !isdefined(@__MODULE__, Symbol("@warn")) # prior to 0.7.0-DEV.2988
    macro warn(args...)
        :(warn($args...))
    end
end

function _indmax(a::AbstractMatrix)
    ind = indmax(a)
    if isa(ind, CartesianIndex) # since 0.7.0-DEV.1660
        return (ind[1], ind[2])
    else
        return ind2sub(size(a), ind)
    end
end

if isdefined(Base, :NamedTuple) # 0.7.0-DEV.2738 to 0.7.0-DEV.3226
    kwargs_pairs(kwargs::NamedTuple) = pairs(kwargs)
end
kwargs_pairs(kwargs) = kwargs
