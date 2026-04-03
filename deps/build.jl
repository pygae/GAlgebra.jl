using PyCall

# Version range for compatible galgebra
const GALGEBRA_MIN = v"0.6.0"
const GALGEBRA_MAX = v"0.7.0"  # exclusive
const GALGEBRA_INSTALL = "0.6.0"

const DEBUG = get(ENV, "GALGEBRA_DEBUG", "") != ""

function debug(msg)
    DEBUG && @info "[GAlgebra.jl] $msg"
end

function detect_galgebra()
    try
        ver_str = pyimport("galgebra").__version__
        ver = VersionNumber(ver_str)
        debug("Found galgebra $ver")
        return ver
    catch e
        debug("galgebra not found: $e")
        return nothing
    end
end

function install_galgebra()
    debug("Installing galgebra==$GALGEBRA_INSTALL via pip")
    run(PyCall.python_cmd(`-m pip install galgebra==$(GALGEBRA_INSTALL) --user`))
end

ver = detect_galgebra()

if ver === nothing
    @info "[GAlgebra.jl] galgebra not found, installing galgebra==$GALGEBRA_INSTALL"
    install_galgebra()
elseif ver < GALGEBRA_MIN || ver >= GALGEBRA_MAX
    @warn "[GAlgebra.jl] galgebra $ver is outside supported range [$GALGEBRA_MIN, $GALGEBRA_MAX). Installing galgebra==$GALGEBRA_INSTALL"
    install_galgebra()
else
    @info "[GAlgebra.jl] Using galgebra $ver (supported range: [$GALGEBRA_MIN, $GALGEBRA_MAX))"
end
