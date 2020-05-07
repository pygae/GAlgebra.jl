# GAlgebra.jl

Julia interface to [GAlgebra](https://github.com/pygae/galgebra), a symbolic Geometric Algebra/Calculus package for SymPy.

## Development Status

Very early. But it already works and has many tests.

| **Build Status**                                                                                |
|:-----------------------------------------------------------------------------------------------:|
| [![](https://travis-ci.com/pygae/GAlgebra.jl.svg?branch=master)](https://travis-ci.com/pygae/GAlgebra.jl) [![](https://img.shields.io/codecov/c/github/pygae/GAlgebra.jl.svg)](https://codecov.io/gh/pygae/GAlgebra.jl) |

## Getting Started

GAlgebra.jl itself doesn't depend on [SymPy.jl](https://github.com/JuliaPy/SymPy.jl), but is designed to work with it.

After installing SymPy.jl and GAlgebra.jl (see below for instructions), you may experiment with GAlgebra.jl just like in the Python version of GAlgebra (though [there're some syntax differences between Python and Julia](https://docs.julialang.org/en/v1/manual/noteworthy-differences/#Noteworthy-differences-from-Python-1) like `True`/`true`, `'`/`"` etc.).

For example, you may start with:

```julia
import SymPy: symbols, sympy, Sym
using GAlgebra

# In console, uncomment to enable colored printing with ANSI escape sequences 
# galgebra.printer.Eprint()
# In Jupyter, uncomment to enable LaTeX printing with MathJax
# galgebra.printer.Format()

(x, y, z) = xyz = symbols("x,y,z",real=true)
o3d = galgebra.ga.Ga("e", g=[1, 1, 1], coords=xyz)

u = o3d.mv("u", "vector")
v = o3d.mv("v", "vector")
A = o3d.mv("A", "mv")
R = o3d.mv("R", "spinor")
# Volume element
I = o3d.I()

# Wedge product: ∧ \wedge
v ∧ A
# Hestenes' inner product: ⋅ \cdot
v ⋅ A
# Left contraction: ⨼ \intprod
v ⨼ A
# Right contraction: ⨽ \intprodr
v ⨽ A
# Scalar product: ⊛ \circledast
# A ⊛ B = <A B†>
v ⊛ A
# Commutator product: ⊠ \boxtimes
# A⊠B = (AB-BA)/2
v ⊠ A
# Anti-commutator product: ⊙ \odot
# A⊙B = (AB+BA)/2
v ⊙ A

# Norm: norm(A) = A.norm() := ||A||
norm(v)

# Inverse: postfix ⁻¹ \^-\^1
# (A)⁻¹ = A^-1 = inv(A) = A.inv()
(R)⁻¹
R^-1
inv(R)

# Reversion: ~A = rev(A) = A.rev()
# A^† is usually used in literature
~A
rev(A)

# Dual: postfix '
# orthogonal complement, Λ^p -> Λ^(n-p)
# note: Ga.dual_mode_value is default to "I+", so A' = A * I
# change Ga.dual_mode_value to get a different definition
A'
dual(A)

# Grade involution: postfix ˣ \^x
# (A)ˣ = A[:*] = involute(A) := A+ - A- = A.even() - A.odd()
# A^* is usually used in literature
(A)ˣ
involute(A)

# Clifford conjugate: postfix ǂ \doublepipe
# (A)ǂ = conj(A) := ((A)^*)^†
(A)ǂ
conj(A)

# Projection: proj(B, A) = A.project_in_blade(B)
proj(u, v)

# Reflection: refl(B, A) = A.reflect_in_blade(B)
refl(u, v)

# Rotation: rot(itheta, A) = A.rotate_multivector(itheta)
# rotate the multivector A by the 2-blade itheta
rot(u ∧ v, A)

# Natural base exponential of x: e^x
exp(u ∧ v)

# Grade-i part: A[i] = A.grade(i) := <A>_i
A[2]

# Scalar (grade-0) part: scalar(A) = A.scalar() := <A> = <A>_0
# note: it returns a SymPy expression unlike A[0] which returns a Mv object
scalar(A)

# Even-grade part: A[:+] = (A)₊ = even(A) = A.even() := A+
A[:+]
even(A)

# Odd-grade part: A[:-] = (A)₋ = odd(A) = A.odd() := A-
A[:-]
odd(A)
```

Note: enter unicode symbols like `∧` with corresponding LaTeX commands like `\wedge` by [Tab completion](https://pkg.julialang.org/docs/julia/THl1k/1.1.0/manual/unicode-input.html) which are provided in the comments.

So far only `galgebra.ga.Ga` and `galgebra.mv.Mv` have been verified to work in Julia, see [tests](https://github.com/pygae/GAlgebra.jl/tree/master/test/runtests.jl). The tests verified many identities in Linear Algebra and Geometric Algebra.

See [examples of GAlgebra](https://github.com/pygae/galgebra/tree/master/examples) for more examples in Python. With some small changes, they can be ported to Julia.

## Installation

The package can be installed with the Julia package manager. From the Julia REPL, type `]` to enter the **Pkg REPL mode** (you'll see a prompt like `(v1.1) pkg>`) and run:

```
dev https://github.com/pygae/GAlgebra.jl.git
```

The installation process will take a while, because it will install [PyCall.jl](https://github.com/JuliaPy/PyCall.jl) and its dependencies and it will also use `pip` to install [GAlgebra](https://github.com/pygae/galgebra) if you haven't. If you wish to use a version other than the latest released version of GAlgebra that will be installed by [deps/build.jl](deps/build.jl), you can install that via pip before/after the installation of `GAlgebra.jl`, e.g. `pip install -e <your local path to GAlgebra>` or `pip install https://github.com/pygae/galgebra/archive/master.zip`.

Now you may run GAlgebra.jl tests in Julia **Pkg REPL mode**: 

```
test GAlgebra
```

At the first time it will take a while, because SymPy.jl is specified as a test dependency so it'll be installed.

Then you'll see something like:

```
   Testing GAlgebra
 Resolving package versions...
Test Summary: | Pass  Broken  Total
GAlgebra.jl   | 1289       1   1290
   Testing GAlgebra tests passed
```

Hint: To get back to the Julia REPL please press backspace, see [Pkg doc](https://docs.julialang.org/en/v1/stdlib/Pkg/index.html) to learn more.
