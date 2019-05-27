module GAlgebra

using PyCall

import Base: convert, getproperty, setproperty!, propertynames

if isdefined(Base, :hasproperty) # Julia 1.2
    import Base: hasproperty
end

import Base: show
import Base: @pure, convert
import Base: +, -, *, /, ^, |, %, ==, !=, <, >, <<, >>, abs, inv, ~, adjoint, getindex

export galgebra
export Mv
# export \cdot, \wedge, \intprod, \intprodr, \dottimes, \timesbar, \circledast
export ⋅, ∧, ⨼, ⨽, ⨰, ⨱, ⊛
# Operator precedence: they have the same precedence, unlike in math
# julia> for op ∈ [:⋅ :∧ :⨼ :⨽ :⨰ :⨱ :⊛]; println(String(op), "  ", Base.operator_precedence(op)) end
# ⋅  13
# ∧  13
# ⨼  13
# ⨽  13
# ⨰  13
# ⨱  13
# ⊛  13

# experimental export \bar\times
# export ×̄

const galgebra = PyCall.PyNULL()
const metric = PyCall.PyNULL()
const ga = PyCall.PyNULL()
const mv = PyCall.PyNULL()
const lt = PyCall.PyNULL()
const printer = PyCall.PyNULL()

mutable struct Mv
    o::PyCall.PyObject
end

convert(::Type{Mv}, o::PyCall.PyObject) = Mv(o)
PyCall.PyObject(o::Mv) = PyCall.PyObject(o.o)

macro define_op(type, op, method)
    @eval begin
        @pure $op(x::$type, y::$type) = x.$method(y)
    end
end

macro define_lop(type, rtype, op, lmethod)
    @eval begin
        @pure $op(x::$type, y::$rtype) = x.$lmethod(y)
    end
end
                
macro define_rop(type, ltype, op, rmethod)
    @eval begin
        @pure $op(x::$ltype, y::$type) = y.$rmethod(x)
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
# Anti-comutator product: \dottimes  A⨰B = (AB+BA)/2
@define_op(Mv, ⨰, __lshift__)
@define_op(Mv, <<, __lshift__)
# # experimental symbol for anti-comutator product: \bar\times A×̄B = (AB+BA)/2
# @define_op(Mv, ×̄, __lshift__)
# Comutator product: \timesbar  A⨱B = (AB-BA)/2
@define_op(Mv, ⨱, __rshift__)
@define_op(Mv, >>, __rshift__)

# Scalar product: \circledast A ⊛ B = <A B†>
@pure ⊛(x::Mv, y::Mv) = (x * ~y).scalar()
@pure %(x::Mv, y::Mv) = x ⊛ y

@pure -(x::Mv) = x.__neg__()

# Norm: abs(A) = |A| = A.norm()
@pure abs(x::Mv) = x.norm()

# Inverse: A^-1 = A.inv()
@pure inv(x::Mv) = x.inv()

# Reverse: ~A = A† = A.rev()
@pure ~(x::Mv) = x.rev()

# Dual: A' = A*
@pure adjoint(x::Mv) = x.dual()

# Grade: A[i] = <A>_i = A.grade(i) = grade-i part of A
@pure getindex(x::Mv, i::Integer) = x.grade(i)

# :⁻¹ => :inv, :⁽²⁾ => :square, :⁽³⁾ => :cube
# :ᵀ => :transpose, :ᴴ => :ctranspose

macro define_postfix_op(super_script, func)
    @eval begin
        struct $super_script end
        export $super_script
        Base.:(*)(x::Mv,::typeof($super_script)) = $func(x)
    end
end

@define_postfix_op(⁻¹, inv)

@define_lop(Mv, Number, +, __add__)
@define_rop(Mv, Number, +, __radd__)
@define_lop(Mv, Number, -, __sub__)
@define_rop(Mv, Number, -, __rsub__)
@define_lop(Mv, Number, *, __mul__)
@define_rop(Mv, Number, *, __rmul__)
@define_lop(Mv, Number, /, __div__)
@define_rop(Mv, Number, /, __rdiv__)
@define_lop(Mv, Number, ==, __eq__)
@define_lop(Mv, Number, !=, __ne__)

@pure function ^(x::Mv, y::Integer)
    if y < 0
        x.__pow__(abs(y)).inv()
    elseif y == 0
        1
    else
        x.__pow__(y)
    end
end

macro define_show(type)
    @eval begin
        Base.show(io::IO, x::$type) = print(io, pystr(x.o))
        Base.show(io::IO, ::MIME"text/plain", x::$type) = print(io, pystr(x.o))
        Base.show(io::IO, ::MIME"text/latex", x::$type) = print(io, "\\begin{align*}" * galgebra.printer.latex(x.o) * "\\end{align*}")
    end
end

@define_show(Mv)

macro delegate_properties(type, obj_field)
    @eval begin
        function getproperty(o::$type, s::AbstractString)
            if s == String($obj_field)
                return getfield(o, $obj_field)
            else
                return getproperty(getfield(o, $obj_field), s)
            end
        end
        
        getproperty(o::$type, s::Symbol) = getproperty(o, String(s))
        
        propertynames(o::$type) = map(x->Symbol(first(x)),
                                        pycall(inspect."getmembers", PyObject, getfield(o, $obj_field)))
        
        # avoiding method ambiguity
        setproperty!(o::$type, s::Symbol, v) = _setproperty!(o,s,v)
        setproperty!(o::$type, s::AbstractString, v) = _setproperty!(o,s,v)
        
        function _setproperty!(o::$type, s::Union{Symbol,AbstractString}, v)
            obj = getfield(o, $obj_field)
            setproperty!(obj, s, v)
            o
        end

        hasproperty(o::$type, s::Symbol) = hasproperty(getfield(o, $obj_field), s)
        hasproperty(o::$type, s::AbstractString) = hasproperty(getfield(o, $obj_field), s)
    end
end

@delegate_properties(Mv, :o)

# Expose Python docstrings to the Julia doc system
Docs.getdoc(x::Mv) = Text(convert(String, x."__doc__"))
Docs.Binding(x::Mv, s::Symbol) = getproperty(x, s)

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
