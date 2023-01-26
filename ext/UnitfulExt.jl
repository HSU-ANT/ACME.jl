# Copyright 2023 Martin Holters
# See accompanying license file.

module UnitfulExt

import ACME
using Unitful: Unitful, @u_str, NoUnits, Quantity, Units, uconvert

remove_unit(unit::Units, q::Number) = NoUnits(uconvert(unit, q) / unit)

ACME.resistor(r::Quantity) = ACME.resistor(remove_unit(u"Ω", r))

ACME.potentiometer(r::Quantity, pos) = ACME.potentiometer(remove_unit(u"Ω", r), pos)
ACME.potentiometer(r::Quantity) = ACME.potentiometer(remove_unit(u"Ω", r))

ACME.capacitor(c::Quantity) = ACME.capacitor(remove_unit(u"F", c))

ACME.inductor(l::Quantity) = ACME.inductor(remove_unit(u"H", l))

function ACME.transformer(
    l1::Quantity, l2::Quantity;
    coupling_coefficient=1,
    mutual_coupling::Quantity=coupling_coefficient*sqrt(l1*l2)
)
    return ACME.transformer(
        remove_unit(u"H", l1), remove_unit(u"H", l2);
        mutual_coupling = remove_unit(u"H", mutual_coupling)
    )
end

function ACME._transformer_ja(ns, α, c, kwargs::NamedTuple{K, <:Tuple{Vararg{Union{Real,Quantity}}}}) where {K}
    units = map(
        (k) -> begin
            if k === :D
                return u"m"
            elseif k === :A
                return u"m^2"
            elseif k === :a || k === :k || k === :Ms
                return u"A/m"
            else
                throw(ArgumentError("transformer: got unsupported keyword argument \"$(k)\""))
            end
        end,
        K
    )
    return ACME._transformer_ja(
        ns, α, c,
        NamedTuple{K}(map((unit, value) -> remove_unit(unit, value), units, values(kwargs))),
    )
end

ACME.voltagesource(v::Quantity; rs::Quantity=0u"Ω") =
    ACME.voltagesource(remove_unit(u"V", v); rs=remove_unit(u"Ω", rs))
ACME._voltagesource(rs::Quantity) = ACME._voltagesource(remove_unit(u"Ω", rs))

ACME.currentsource(i::Quantity; gp::Quantity=0u"S") =
    ACME.currentsource(remove_unit(u"A", i); gp=remove_unit(u"S", gp))
ACME._currentsource(gp::Quantity) = ACME._currentsource(remove_unit(u"S", gp))

ACME._voltageprobe(gp::Quantity) = ACME._voltageprobe(remove_unit(u"S", gp))

ACME._currentprobe(rs::Quantity) = ACME._currentprobe(remove_unit(u"Ω", rs))

ACME._diode(is::Quantity, η::Real) = ACME._diode(remove_unit(u"A", is), η)

function ACME._bjt(typ, kwargs::NamedTuple{K, <:Tuple{Vararg{Union{Real,Quantity}}}}) where {K}
    units = map(
        (k) -> begin
            if k === :is || k === :isc || k === :ise || k === :ilc || k === :ile || k === :ikf || k === :ikr
                return u"A"
            elseif k === :vaf || k === :var
                return u"V"
            elseif k === :re || k === :rc || k === :rb
                return u"Ω"
            elseif k === :η || k === :ηc || k === :ηe || k === :βf || k === :βr || k === :ηcl || k === :ηel
                return NoUnits
            else
                throw(ArgumentError("bjt: got unsupported keyword argument \"$(k)\""))
            end
        end,
        K
    )
    return ACME._bjt(
        typ,
        NamedTuple{K}(map((unit, value) -> remove_unit(unit, value), units, values(kwargs))),
    )
end

_mosfet_remove_units(unit::Units, q::Number) = remove_unit(unit, q)
_mosfet_remove_units(unit::Units, q::Tuple{Vararg{Number,N}}) where {N} =
    ntuple(n -> remove_unit(unit / u"V"^(n-1), q[n]), Val(N))

function ACME._mosfet(typ, kwargs::NamedTuple{K, <:Tuple{Vararg{Union{Union{Real,Quantity},Tuple{Vararg{Union{Real,Quantity}}}}}}}) where {K}
    units = map(
        (k) -> begin
            if k === :vt
                return u"V"
            elseif k === :α
                return u"A/V^2"
            elseif k === :λ
                return u"V^-1"
            else
                throw(ArgumentError("bjt: got unsupported keyword argument \"$(k)\""))
            end
        end,
        K
    )
    return ACME._mosfet(
        typ,
        NamedTuple{K}(map((unit, q) -> _mosfet_remove_units(unit, q), units, values(kwargs))),
    )
end

ACME.opamp(::Type{Val{:macak}}, gain::Real, vomin::Quantity, vomax::Quantity) =
    ACME.opamp(Val{:macak}, gain, remove_unit(u"V", vomin), remove_unit(u"V", vomax))

end
