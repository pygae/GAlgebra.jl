module GAlgebra

using PyCall

import Base: @pure 

export galgebra
export Mv

const galgebra = PyCall.PyNULL()
const metric = PyCall.PyNULL()
const ga = PyCall.PyNULL()
const mv = PyCall.PyNULL()
const lt = PyCall.PyNULL()
const printer = PyCall.PyNULL()

if isdefined(Base, :hasproperty) # Julia 1.2
    import Base: hasproperty
end

macro define_show(type)
    @eval begin
        Base.show(io::IO, x::$type) = print(io, pystr(x.o))
        Base.show(io::IO, ::MIME"text/plain", x::$type) = print(io, pystr(x.o))
        Base.show(io::IO, ::MIME"text/latex", x::$type) = print(io, "\\begin{align*}" * galgebra.printer.latex(x.o) * "\\end{align*}")
    end
end

macro delegate_properties(type, obj_field)
    @eval begin
        Base.convert(::Type{$type}, o::PyCall.PyObject) = $type(o)
        PyCall.PyObject(o::$type) = PyCall.PyObject(getfield(o, $obj_field))

        function Base.getproperty(o::$type, s::AbstractString)
            if s == String($obj_field)
                return getfield(o, $obj_field)
            else
                return getproperty(getfield(o, $obj_field), s)
            end
        end
        
        Base.getproperty(o::$type, s::Symbol) = getproperty(o, String(s))
        
        Base.propertynames(o::$type) = map(x->Symbol(first(x)),
                                        pycall(inspect."getmembers", PyObject, getfield(o, $obj_field)))
        
        # avoiding method ambiguity
        Base.setproperty!(o::$type, s::Symbol, v) = _setproperty!(o,s,v)
        Base.setproperty!(o::$type, s::AbstractString, v) = _setproperty!(o,s,v)
        
        function _setproperty!(o::$type, s::Union{Symbol,AbstractString}, v)
            obj = getfield(o, $obj_field)
            setproperty!(obj, s, v)
            o
        end

        hasproperty(o::$type, s::Symbol) = hasproperty(getfield(o, $obj_field), s)
        hasproperty(o::$type, s::AbstractString) = hasproperty(getfield(o, $obj_field), s)
    end
end

macro delegate_doc(type)
    @eval begin
        # Expose Python docstrings to the Julia doc system
        Docs.getdoc(x::$type) = Text(convert(String, x."__doc__"))
        Docs.Binding(x::$type, s::Symbol) = getproperty(x, s)
    end
end

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

macro define_unary_op(type, op, method)
    @eval begin
        @pure $op(x::$type) = x.$method()
    end
end

# :⁻¹ => :inv, :⁽²⁾ => :square, :⁽³⁾ => :cube
# :ᵀ => :transpose, :ᴴ => :ctranspose

macro define_postfix_symbol(super_script)
    @eval begin
        const $super_script = () -> String(:super_script)
        export $super_script
    end
end

# \^-\^1
@define_postfix_symbol(⁻¹)
# \^T
@define_postfix_symbol(ᵀ)
# \doublepipe
@define_postfix_symbol(ǂ)
# \^-
@define_postfix_symbol(⁻)

macro define_postfix_op(type, super_script, func)
    @eval begin
        Base.:(*)(x::$type,::typeof($super_script)) = $func(x)
        Base.:(^)(x::$type,::typeof($super_script)) = $func(x)
    end
end

mutable struct Mv
    o::PyCall.PyObject
end

@define_show(Mv)
@delegate_properties(Mv, :o)
@delegate_doc(Mv)

import Base: +, -, *, /, ^, |, %, ==, !=, <, >, <<, >>, abs, inv, ~, adjoint, getindex
# export \cdot, \wedge, \intprod, \intprodr, \dottimes, \timesbar, \circledast, \times
export ⋅, ∧, ⨼, ⨽, ⨰, ⨱, ⊛, ×
# Operator precedence: they have the same precedence, unlike in math
# julia> for op ∈ [:* :⋅ :∧ :⨼ :⨽ :⨰ :⨱ :⊛ :×]; println(String(op), "  ", Base.operator_precedence(op)) end
# *  13
# ⋅  13
# ∧  13
# ⨼  13
# ⨽  13
# ⨰  13
# ⨱  13
# ⊛  13
# ×  13

# experimental export \bar\times
# export ×̄

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
# Anti-comutator product: \dottimes
# A⨰B = (AB+BA)/2
@define_op(Mv, ⨰, __lshift__)
@define_op(Mv, <<, __lshift__)

# Cross product for vectors in 3D
@pure ×(x::Mv, y::Mv) = mv.cross(x, y)

# # experimental symbol for anti-comutator product: \bar\times
# A×̄B = (AB+BA)/2
# @define_op(Mv, ×̄, __lshift__)

# Comutator product: \timesbar
# A⨱B = (AB-BA)/2
@define_op(Mv, ⨱, __rshift__)
@define_op(Mv, >>, __rshift__)

# Scalar product: \circledast
# A ⊛ B = <A B†>
@pure ⊛(x::Mv, y::Mv) = (x * ~y).scalar()
@pure %(x::Mv, y::Mv) = x ⊛ y

@define_unary_op(Mv, -, __neg__)

# Norm: abs(A) = |A| = A.norm()
@define_unary_op(Mv, abs, norm)

# Inverse: \^-\^1
# (A)⁻¹ = A^-1 = A.inv()
@define_unary_op(Mv, inv, inv)
@define_postfix_op(Mv, ⁻¹, inv)

# Reversion: ~A = (A)ᵀ = A.rev()
# A^† is usually used in literature, but \dagger is reserved by Julia
@define_unary_op(Mv, ~, rev)
@define_unary_op(Mv, rev, rev)

# @deprecated
# @define_postfix_op(Mv, ᵀ, rev)

# Dual: A' = A * I
# note: Ga.dual_mode_value is default to "I+"
# change Ga.dual_mode_value to get a different definition
# A^⊥ (\bot) is sometimes used in literature
@define_unary_op(Mv, adjoint, dual)

export rev, involute, proj, refl, rot, exp_with_hint

# Grade involution: \^-
# (A)⁻ = A+ - A-
# A^* is usually used in literature but * and ∗(\ast) both parsed as binary operator in Julia
@pure involute(x::Mv) = x.even() - x.odd()
@define_postfix_op(Mv, ⁻, involute)

# Clifford conjugate: \doublepipe
# (A)ǂ = ((A)^*)^†
# A^‡ is usually used in literature but \ddagger is reserved by Julia
@pure Base.conj(x::Mv) = involute(x).rev()
@define_postfix_op(Mv, ǂ, conj)

# Grade-i part: A[i] = <A>_i = A.grade(i)
@pure getindex(x::Mv, i::Integer) = x.grade(i)

# Projection: proj(B, A) = A.project_in_blade(B)
@pure proj(y::Mv, x::Mv) = mv.proj(y, x)

# Reflection: refl(B, A) = A.reflect_in_blade(B)
@pure refl(y::Mv, x::Mv) = mv.refl(y, x)

# Rotate by the 2-blade itheta the multivector A
# rot(itheta, A, hint="-") = A.rotate_multivector(itheta, hint)
@pure rot(itheta::Mv, A::Mv, hint::AbstractString="-") = mv.rot(itheta, A, hint)

# Natural base exponential of x: e^x
@pure Base.exp(x::Mv) = exp_with_hint(x, "-")
@pure exp_with_hint(x::Mv, hint::AbstractString="-") = mv.exp(x, hint)

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
