# Copyright 2018 Martin Holters
# See accompanying license file.

import Compat
using Compat: @warn, BitSet, argmax, argmin, copyto!, dropdims, findall, rmul!, undef
using Compat.Markdown: @doc_str

if isdefined(Base, :NamedTuple) # 0.7.0-DEV.2738 to 0.7.0-DEV.3226
    kwargs_pairs(kwargs::NamedTuple) = pairs(kwargs)
end
kwargs_pairs(kwargs) = kwargs

@static if !isdefined(Compat.SparseArrays, :blockdiag) # prior to 0.7.0-DEV.4498
    const blockdiag = Compat.SparseArrays.blkdiag
else
    using Compat.SparseArrays: blockdiag
end
