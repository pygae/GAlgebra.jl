# operators imported and overridden from Base
import Base: +, -, *, /, ^, |, %, ==, !=, <, >, <<, >>, ~

# functions imported and overridden from Base
# convention: prefix `Base.` should always be present when overriding
# so comment the following will still work
import Base: show,convert, getproperty, propertynames, setproperty!
if isdefined(Base, :hasproperty) # Julia 1.2
    import Base: hasproperty
end
import Base: abs, inv, adjoint, exp, getindex

export Mv

export norm, rev, dual, involute, proj, refl, rot, exp_with_hint, scalar, even, odd
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

# \^-\^1
@define_postfix_symbol(⁻¹)
# # \^T
# @define_postfix_symbol(ᵀ)
# \^x
@define_postfix_symbol(ˣ)
# \doublepipe
@define_postfix_symbol(ǂ)
# \_+
@define_postfix_symbol(₊)
# \_-
@define_postfix_symbol(₋)
# \bot
# @define_postfix_symbol(⊥)

mutable struct Mv
    o::PyCall.PyObject
end

@define_show(Mv)
@delegate_properties(Mv, :o)
@delegate_doc(Mv)

@define_op(Mv, +, __add__)
@define_op(Mv, -, __sub__)
# Geometric product: *
@define_op(Mv, *, __mul__)
@define_op(Mv, /, __div__)
@define_op(Mv, ==, __eq__)
@define_op(Mv, !=, __ne__)

# Wedge product: \wedge
@define_op(Mv, ∧, __xor__)
# Hestenes' inner product: ⋅ \cdot
@define_op(Mv, ⋅, __or__)
@define_op(Mv, |, __or__)
# Left contraction: \intprod
@define_op(Mv, ⨼, __lt__)
@define_op(Mv, <, __lt__)
# Right contraction: \intprodr
@define_op(Mv, ⨽, __gt__)
@define_op(Mv, >, __gt__)

# Commutator product: ⨱ \timesbar
# A⨱B = (AB-BA)/2
@define_op(Mv, ⨱, __rshift__)
@define_op(Mv, >>, __rshift__)

# Anti-commutator product: ⨰ \dottimes
# A⨰B = (AB+BA)/2
@define_op(Mv, ⨰, __lshift__)
@define_op(Mv, <<, __lshift__)

# # experimental symbol for anti-comutator product: \bar\times
# A×̄B = (AB+BA)/2
# @define_op(Mv, ×̄, __lshift__)

# Cross product for vectors in 3D
@pure ×(x::Mv, y::Mv) = mv.cross(x, y)

# Scalar product: \circledast
# A ⊛ B = <A B†>
@pure ⊛(x::Mv, y::Mv) = (x * ~y).scalar()
@pure %(x::Mv, y::Mv) = x ⊛ y

@define_unary_op(Mv, -, __neg__)

# Norm: abs(A) = norm(A) = A.norm() := ||A||
@define_unary_op(Mv, Base.abs, norm)
@define_unary_op(Mv, norm, norm)

# Inverse: \^-\^1
# (A)⁻¹ = A^-1 = A.inv()
@define_unary_op(Mv, Base.inv, inv)
@define_postfix_op(Mv, ⁻¹, Base.inv)

# Reversion: ~A = A[:~] = A.rev()
# A^† is usually used in literature, but \dagger is reserved by Julia
@define_unary_op(Mv, ~, rev)
@define_unary_op(Mv, rev, rev)

# @deprecated
# @define_postfix_op(Mv, ᵀ, rev)

# Dual: A' = A * I
# note: Ga.dual_mode_value is default to "I+"
# change Ga.dual_mode_value to get a different definition
# A^⊥ (\bot) is sometimes used in literature
@define_unary_op(Mv, Base.adjoint, dual)
@define_unary_op(Mv, dual, dual)
# @define_postfix_op(Mv, ⊥, dual)

# Grade involution: postfix ˣ \^x
# (A)ˣ = A[:*] = involute(A) := A+ - A- = A.even() - A.odd()
# A^* is usually used in literature
@pure involute(x::Mv) = x.even() - x.odd()
@define_postfix_op(Mv, ˣ, involute)

# Clifford conjugate: \doublepipe
# (A)ǂ = A[:ǂ] := ((A)^*)^†
# A^‡ is usually used in literature but \ddagger is reserved by Julia
@pure Base.conj(x::Mv) = involute(x).rev()
@define_postfix_op(Mv, ǂ, Base.conj)

# Projection: proj(B, A) = A.project_in_blade(B)
@pure proj(y::Mv, x::Mv) = mv.proj(y, x)

# Reflection: refl(B, A) = A.reflect_in_blade(B)
@pure refl(y::Mv, x::Mv) = mv.refl(y, x)

# Rotation: rot(itheta, A) = A.rotate_multivector(itheta)
# rotate the multivector A by the 2-blade itheta
@pure rot(itheta::Mv, A::Mv, hint::AbstractString="-") = mv.rot(itheta, A, hint)

# Natural base exponential of x: e^x
@pure Base.exp(x::Mv) = exp_with_hint(x, "-")
@pure exp_with_hint(x::Mv, hint::AbstractString="-") = mv.exp(x, hint)

# Grade-i part: A[i] = <A>_i = A.grade(i)
@pure Base.getindex(x::Mv, i::Integer) = x.grade(i)
@pure function Base.getindex(x::Mv, sym::Symbol)
    if sym == :+
        even(x)
    elseif sym == :-
        odd(x)
    elseif sym == :~
        rev(x)
    elseif sym == :⁻¹
        inv(x)
    elseif sym == :*
        involute(x)
    elseif sym == :ǂ
        conj(x)
    else
        throw(DomainError(x, "argument can only be one of :+, :-, :~, :⁻¹, :*, :ǂ"))
    end
end

# Scalar (grade-0) part: scalar(A) = A.scalar() := <A> = <A>_0
# note: it returns a SymPy expression unlike A[0] which returns a Mv object
@define_unary_op(Mv, scalar, scalar)

# Even-grade part: A[:+] = (A)₊ = even(A) = A.even() := A+
@define_unary_op(Mv, even, even)
@define_postfix_op(Mv, ₊, even)

# Odd-grade part: A[:-] = (A)₋ = odd(A) = A.odd() := A-
@define_unary_op(Mv, odd, odd)
@define_postfix_op(Mv, ₋, odd)

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
