# Copyright 2019 Martin Holters
# See accompanying license file.

function wrap_nleq_expr(nn, nq, expr)
    Base.depwarn("nonlinear_eq should be given as a function, not an expression", :Element)

    res_symbs = [gensym("res_$n") for n in 1:nn]
    J_symbs = [gensym("J_$(m)_$n") for m in 1:nn, n in 1:nq]

    function rewrite_refs(expr::Expr)
        if expr.head == :ref
            if expr.args[1] === :res
                return res_symbs[expr.args[2]]
            elseif expr.args[1] === :J
                return J_symbs[expr.args[2], expr.args[3]]
            end
            return expr
        else
            return Expr(expr.head, rewrite_refs.(expr.args)...)
        end
    end

    rewrite_refs(x::Any) = x

    expr = rewrite_refs(expr)

    return eval(quote
        @inline function (q)
            $(nn > 0 || nq > 0 ? Expr(:local, res_symbs..., J_symbs...) : nothing)
            $(expr)
            res = SVector{$nn}($(res_symbs...))
            J = SMatrix{$nn,$nq}($(J_symbs...))
            return (res, J)
        end
    end)
end

function ports_from_old_pins(pins)
    ports = [pins[2i-1] => pins[2i] for i in 1:length(pins)÷2]
    Base.depwarn("`Element(..., pins=$(repr(pins)))` is depreated, use `Element(..., ports=$(repr(ports, context=:typeinfo=>typeof(ports))))` instead", :Element)
    return ports
end

Base.@deprecate(composite_element(circ::Circuit, pins::Vector{<:Pair}),
    composite_element(circ,  pinmap=Dict(pins...)))
