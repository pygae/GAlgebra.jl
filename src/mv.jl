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
# export \cdot, \wedge, \intprod, \intprodr, \odot, \boxtimes, \circledast, \times
export +, -, *, /, ^, |, %, ==, !=, <, >, <<, >>, ~, ⋅, ∧, ⨼, ⨽, ⊙, ⊠, ⊛, ×
# Operator precedence: they have the same precedence, unlike in math
# julia> for op ∈ [:* :⋅ :∧ :⨼ :⨽ :⊙ :⊠ :⊛ :×]; println(String(op), "  ", Base.operator_precedence(op)) end
# *  13
# ⋅  13
# ∧  13
# ⨼  13
# ⨽  13
# ⊙  13
# ⊠  13
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

@doc raw"""
A wrapper class for `galgebra.mv.Mv`:
    
- all methods of `galgebra.mv.Mv` are delegated and can be called like in Python.
- enhanced with operator overriding and extra methods.
- automatically supports pretty printing like in Python.
"""
mutable struct Mv
    o::PyCall.PyObject
end

@define_show(Mv)
@delegate_properties(Mv, :o)
@delegate_doc(Mv)

"Addition."
@define_op(Mv, +, __add__)

"Subtraction."
@define_op(Mv, -, __sub__)

@doc raw"""
Geometric product.

``A * B \equiv A B``.
"""
@define_op(Mv, *, __mul__)

@doc raw"""
Division.

``A / B \equiv A B^{-1}``. Only valid when ``B`` has inverse.
"""
@define_op(Mv, /, __div__)

@doc raw"""
Comparisons of equality.

``A = B \equiv \mathrm{simplify}(A - B) = 0``
"""
@define_op(Mv, ==, __eq__)

@doc raw"""
Comparisons of inequality.

Hint: Type ≠ with `\neq`. Alternatively, use `!=`.

``A \neq B \equiv \mathrm{simplify}(A - B) \neq 0``
"""
@define_op(Mv, ≠, __ne__)
@define_op(Mv, !=, __ne__)

@doc raw"""
Wedge product.

Hint: type ∧ with `\wedge`.
"""
@define_op(Mv, ∧, __xor__)

@doc raw"""
Hestenes' inner product.

Hint: type ⋅ with `\cdot`. Alternatively, use `|`.

``\boldsymbol{A}_{r} \cdot \boldsymbol{B}_{s} \equiv \left\{\begin{array}{lr}{r \text { and } s \neq 0 :} & {\left\langle\boldsymbol{A}_{r} \boldsymbol{B}_{s}\right\rangle_{|r-s|}} \\ {r \text { or } s=0 :} & {0}\end{array}\right.``
"""
@define_op(Mv, ⋅, __or__)
@define_op(Mv, |, __or__)

@doc raw"""Left contraction, i.e. "contraction onto". 

Hint: type ⨼ with `\intprod`. Alternatively, use `<`.

``A \rfloor B \equiv \sum\limits_{r, s}\left\langle\langle A\rangle_{r}\langle B\rangle_{s}\right\rangle_{s-r}``

In literature the notation is usually ``A \rfloor B``, but `\rfloor` is reserved by Julia.
"""
@define_op(Mv, ⨼, __lt__)
@define_op(Mv, <, __lt__)

@doc raw"""Right contraction, i.e. "contraction by".

Hint: type ⨽ with `\intprodr`. Alternatively, use `>`.

``A \lfloor B \equiv \sum\limits_{r, s}\left\langle\langle A\rangle_{r}\langle B\rangle_{s}\right\rangle_{r-s}``

In literature the notation is usually ``A \lfloor B``, but `\lfloor` is reserved by Julia.
"""
@define_op(Mv, ⨽, __gt__)
@define_op(Mv, >, __gt__)


@doc raw"""Commutator product.

Hint: type ⊠ with `\boxtimes`. Alternatively, use `>>`.

``A \underline{\times} B \equiv \dfrac{1}{2}(AB-BA)``.
"""
@define_op(Mv, ⊠, __rshift__)
@define_op(Mv, >>, __rshift__)

@doc raw"""Anti-commutator product.

Hint: type ⊙ with `\odot`. Alternatively, use `<<`.

``A \bar{\times} B \equiv \dfrac{1}{2}(AB+BA)``.
"""
@define_op(Mv, ⊙, __lshift__)
@define_op(Mv, <<, __lshift__)

# # experimental symbol for anti-comutator product: \bar\times
# A×̄B = (AB+BA)/2
# @define_op(Mv, ×̄, __lshift__)

@doc raw"""
Cross product for vectors in 3D.
"""
@define_op_with_impl(Mv, ×, mv.cross(x, y))

@doc raw"""Scalar product.

Hint: type ⊛ with `\circledast`. Alternatively, use `%`.

``A \circledast B \equiv \langle A B^{\dagger} \rangle``.

In literature the notation is usually ``\ast`` , but it's visually indistinguishable from `*`.
"""
@define_op_with_impl(Mv, ⊛, (x * ~y).scalar())
@define_op_with_impl(Mv, %, x ⊛ y)

"Unary negation."
@define_unary_op(Mv, -, __neg__)

@doc raw"""
Norm.

`norm(A)` = `A.norm()` ``\equiv \left\lVert A \right\rVert \equiv \sqrt{A \tilde{A}}``

Alternatively:

- `A.norm(hint="+")` ``\equiv \sqrt{A \tilde{A}}``
- `A.norm(hint="-")` ``\equiv \sqrt{- A \tilde{A}}``
- `A.norm(hint="0")` ``\equiv \sqrt{\left| A \tilde{A} \right|}``

Only valid when the result is a scalar.
"""
@define_unary_op(Mv, norm, norm)
@define_unary_op(Mv, Base.abs, norm)

@doc raw"""
Inverse.

`(A)⁻¹ = A^-1 = inv(A) = A.inv()` ``\equiv A^{-1}``

Hint: type ⁻¹ with `\^-\^1`.
"""
@define_unary_op(Mv, Base.inv, inv)
@define_postfix_op(Mv, ⁻¹, Base.inv)

@doc raw"""
Reversion.

`~A = A[:~] = rev(A) = A.rev()` ``\equiv \tilde{A} \equiv A^{\dagger}``

In literature the notation is usually ``\tilde{A}`` or ``A^{\dagger}``, the former is illegal syntax and `\dagger` in the latter is is reserved by Julia.
"""
@define_unary_op(Mv, ~, rev)
@define_unary_op(Mv, rev, rev)

# @deprecated
# @define_postfix_op(Mv, ᵀ, rev)

@doc raw"""
Dual, i.e. orthogonal complement, ``\Lambda^p \to \Lambda^{n-p}``.

`A'` ``\equiv A^{\bot} \equiv A I``

Note: call `Ga.dual_mode(mode)` to globally specify a different dual mode (`I+` is the default):

| dual_mode | ``A^{\bot}`` |
|-----------|--------------|
|   `+I`    |     ``IA``   |
|   `-I`    |    ``-IA``   |
|   `I+`    |     ``AI``   |
|   `I-`    |    ``-AI``   |
|  `+Iinv`  |  ``I^{-1}A`` |
|  `-Iinv`  | ``-I^{-1}A`` |
|  `Iinv+`  |  ``AI^{-1}`` |
|  `Iinv-`  | ``-AI^{-1}`` |
"""
@define_unary_op(Mv, Base.adjoint, dual)
@define_unary_op(Mv, dual, dual)
# @define_postfix_op(Mv, ⊥, dual)

@doc raw"""
Grade involution.

`(A)ˣ = A[:*] = involute(A)` ``\equiv A_+ - A_- \equiv`` `A.even() - A.odd()`

Hint: type ˣ with `\^x`.

In literature the notation is usually ``A^{*}``.
"""
@define_unary_op_with_impl(Mv, involute, x.even() - x.odd())
@define_postfix_op(Mv, ˣ, involute)

@doc raw"""
Clifford conjugate.

`(A)ǂ = A[:ǂ]` ``\equiv A^{*\dagger}``

Hint: type ǂ with `\doublepipe`.

In literature the notation is usually ``A^{\ddagger}``, but `\ddagger` is reserved by Julia.
"""
@define_unary_op_with_impl(Mv, Base.conj, involute(x).rev())
@define_postfix_op(Mv, ǂ, Base.conj)

@doc raw"""
Projection.

`proj(B, A)` `` \equiv P_{B}(A) \equiv`` `A.project_in_blade(B)`

Only valid if B is a blade.
"""
@define_op_with_impl(Mv, proj, mv.proj(x, y))

@doc raw"""
Reflection.

`refl(B, A)` `` \equiv \mathrm{Refl}_{B}(A) \equiv`` `A.reflect_in_blade(B)`

Only valid if B is a blade.
"""
@define_op_with_impl(Mv, refl, mv.refl(x, y))

@doc raw"""
Rotation.

Rotate the multivector `A` by the 2-blade `itheta`.

`rot(itheta, A)` ``\equiv A e^{I \theta} \equiv`` `A.rotate_multivector(itheta)`
"""
@define_op_with_impl(Mv, rot, rot_with_hint(x, y))
@pure rot_with_hint(itheta::Mv, A::Mv, hint::AbstractString="-") = mv.rot(itheta, A, hint)

@doc raw"Natural base exponential of X: ``e^X``"
@define_unary_op_with_impl(Mv, Base.exp, exp_with_hint(x))
@pure exp_with_hint(x::Mv, hint::AbstractString="-") = mv.exp(x, hint)

@doc raw"""
The `i`-th grade part.

`A[i] = A.grade(i)` ``\equiv \langle A B^{\dagger} \rangle_i``
"""
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

@doc raw"""
Scalar (grade-0) part.

`scalar(A) = A.scalar()` ``\equiv \langle A B^{\dagger} \rangle \equiv \langle A B^{\dagger} \rangle_0``

Note: it returns a SymPy expression unlike A[0] which returns a Mv object
"""
@define_unary_op(Mv, scalar, scalar)

@doc raw"""
Even-grade part.

`A[:+] = (A)₊ = even(A) = A.even()` ``\equiv A_+``
"""
@define_unary_op(Mv, even, even)
@define_postfix_op(Mv, ₊, even)

@doc raw"""
Odd-grade part.

`A[:-] = (A)₋ = odd(A) = A.odd()` ``\equiv A_-``
"""
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
