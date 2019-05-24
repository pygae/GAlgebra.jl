module GAlgebra

using PyCall

import Base: show
import Base: convert
import Base: +,-,*,/,^,|,==,!=,<,>,<<,>>

export galgebra
export Mv
# export \cdot, \wedge, \intprod, \intprodr, \timesbar
export ⋅,∧,⨼,⨽,⨱

const galgebra = PyCall.PyNULL()
const metric = PyCall.PyNULL()
const ga = PyCall.PyNULL()
const mv = PyCall.PyNULL()
const lt = PyCall.PyNULL()
const printer = PyCall.PyNULL()

mutable struct Mv
    o::PyCall.PyObject
end

Base.convert(::Type{Mv}, o::PyCall.PyObject) = Mv(o)

macro define_op(type, op, method)
    @eval begin
        $op(x::$type, y::$type) = x.o.$method(y.o)
    end
end

macro define_lop(type, rtype, op, lmethod)
    @eval begin
        $op(x::$type, y::$rtype) = x.o.$lmethod(y)
    end
end
                
macro define_rop(type, ltype, op, rmethod)
    @eval begin
        $op(x::$ltype, y::$type) = y.o.$rmethod(x)
    end
end

@define_op(Mv, +, __add__)
@define_op(Mv, -, __sub__)
# Geometric product: *
@define_op(Mv, *, __mul__)
@define_op(Mv, /, __div__)
@define_op(Mv, ==, __eq__)
@define_op(Mv, !=, __ne__)

# Wedge product: \wedge
@define_op(Mv, ∧, __xor__)
# Hestene's inner product: \cdot
@define_op(Mv, ⋅, __or__)
@define_op(Mv, |, __or__)
# Left contraction: \intprod
@define_op(Mv, ⨼, __lt__)
@define_op(Mv, <, __lt__)
# Right contraction: \intprodr
@define_op(Mv, ⨽, __gt__)
@define_op(Mv, >, __gt__)
# Anti-comutator product: A<<B = (AB+BA)/2
@define_op(Mv, <<, __lshift__)
# Comutator product: \timesbar  A⨱B = (AB-BA)/2
@define_op(Mv, ⨱, __rshift__)
@define_op(Mv, >>, __rshift__)

@define_lop(Mv, Number, +, __add__)
@define_rop(Mv, Number, +, __radd__)
@define_lop(Mv, Number, -, __sub__)
@define_rop(Mv, Number, -, __rsub__)
@define_lop(Mv, Number, *, __mul__)
@define_rop(Mv, Number, *, __rmul__)
@define_lop(Mv, Number, /, __div__)
@define_rop(Mv, Number, /, __rdiv__)
@define_lop(Mv, Number, ^, __pow__)
@define_lop(Mv, Number, ==, __eq__)
@define_lop(Mv, Number, !=, __ne__)

-(x::Mv) = x.o.__neg__()

macro define_show(type)
    @eval begin
        Base.show(io::IO, x::$type) = print(io, pystr(x.o))
        Base.show(io::IO, ::MIME"text/plain", x::$type) = print(io, pystr(x.o))
        Base.show(io::IO, ::MIME"text/latex", x::$type) = print(io, "\\begin{align*}" * galgebra.printer.latex(x.o) * "\\end{align*}")
    end
end

@define_show(Mv)

function __init__()
    copy!(galgebra, PyCall.pyimport_conda("galgebra", "galgebra"))
    copy!(metric, PyCall.pyimport_conda("galgebra.metric", "galgebra"))
    copy!(ga, PyCall.pyimport_conda("galgebra.ga", "galgebra"))
    copy!(mv, PyCall.pyimport_conda("galgebra.mv", "galgebra"))
    copy!(lt, PyCall.pyimport_conda("galgebra.lt", "galgebra"))
    copy!(printer, PyCall.pyimport_conda("galgebra.printer", "galgebra"))
    
    pytype_mapping(galgebra.mv.Mv, Mv)
end

end # module
