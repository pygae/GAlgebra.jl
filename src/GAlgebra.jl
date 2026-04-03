module GAlgebra

using PyCall

import Base: @pure

export galgebra

const galgebra = PyCall.PyNULL()
const metric = PyCall.PyNULL()
const ga = PyCall.PyNULL()
const mv = PyCall.PyNULL()
const lt = PyCall.PyNULL()
const printer = PyCall.PyNULL()

include("macros.jl")
include("ga.jl")
include("mv.jl")

function __init__()
    copy!(galgebra, pyimport("galgebra"))
    copy!(metric, pyimport("galgebra.metric"))
    copy!(ga, pyimport("galgebra.ga"))
    copy!(mv, pyimport("galgebra.mv"))
    copy!(lt, pyimport("galgebra.lt"))
    copy!(printer, pyimport("galgebra.printer"))

    if get(ENV, "GALGEBRA_DEBUG", "") != ""
        ver = galgebra.__version__
        python = PyCall.pyprogramname
        @info "[GAlgebra.jl] galgebra $ver loaded from Python: $python"
    end

    pytype_mapping(galgebra.mv.Mv, Mv)
end

end # module
