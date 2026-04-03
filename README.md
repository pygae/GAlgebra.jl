# GAlgebra.jl

Julia interface to [GAlgebra](https://github.com/pygae/galgebra), a symbolic Geometric Algebra/Calculus package for SymPy.

## Development Status

Very early. But it already works and has many tests.

| **Build Status**                                                                                |
|:-----------------------------------------------------------------------------------------------:|
| ![CI](https://github.com/pygae/GAlgebra.jl/workflows/CI/badge.svg) [![](https://img.shields.io/codecov/c/github/pygae/GAlgebra.jl.svg)](https://codecov.io/gh/pygae/GAlgebra.jl) |

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

# Wedge product: â§ \wedge
v â§ A
# Hestenes' inner product: â \cdot
v â A
# Left contraction: â¨¼ \intprod
v â¨¼ A
# Right contraction: â¨½ \intprodr
v â¨½ A
# Scalar product: â \circledast
# A â B = <A Bâ >
v â A
# Commutator product: â  \boxtimes
# Aâ B = (AB-BA)/2
v â  A
# Anti-commutator product: â \odot
# AâB = (AB+BA)/2
v â A

# Norm: norm(A) = A.norm() := ||A||
norm(v)

# Inverse: postfix â»Â¹ \^-\^1
# (A)â»Â¹ = A^-1 = inv(A) = A.inv()
(R)â»Â¹
R^-1
inv(R)

# Reversion: ~A = rev(A) = A.rev()
# A^â  is usually used in literature
~A
rev(A)

# Dual: postfix '
# orthogonal complement, Î^p -> Î^(n-p)
# note: Ga.dual_mode_value is default to "I+", so A' = A * I
# change Ga.dual_mode_value to get a different definition
A'
dual(A)

# Grade involution: postfix Ë£ \^x
# (A)Ë£ = A[:*] = involute(A) := A+ - A- = A.even() - A.odd()
# A^* is usually used in literature
(A)Ë£
involute(A)

# Clifford conjugate: postfix Ç \doublepipe
# (A)Ç = conj(A) := ((A)^*)^â 
(A)Ç
conj(A)

# Projection: proj(B, A) = A.project_in_blade(B)
proj(u, v)

# Reflection: refl(B, A) = A.reflect_in_blade(B)
refl(u, v)

# Rotation: rot(itheta, A) = A.rotate_multivector(itheta)
# rotate the multivector A by the 2-blade itheta
rot(u â§ v, A)

# Natural base exponential of x: e^x
exp(u â§ v)

# Grade-i part: A[i] = A.grade(i) := <A>_i
A[2]

# Scalar (grade-0) part: scalar(A) = A.scalar() := <A> = <A>_0
# note: it returns a SymPy expression unlike A[0] which returns a Mv object
scalar(A)

# Even-grade part: A[:+] = (A)â = even(A) = A.even() := A+
A[:+]
even(A)

# Odd-grade part: A[:-] = (A)â = odd(A) = A.odd() := A-
A[:-]
odd(A)
```

Note: enter unicode symbols like `â§` with corresponding LaTeX commands like `\wedge` by [Tab completion](https://pkg.julialang.org/docs/julia/THl1k/1.1.0/manual/unicode-input.html) which are provided in the comments.

So far only `galgebra.ga.Ga` and `galgebra.mv.Mv` have been verified to work in Julia, see [tests](https://github.com/pygae/GAlgebra.jl/tree/master/test/runtests.jl). The tests verified many identities in Linear Algebra and Geometric Algebra.

See [examples of GAlgebra](https://github.com/pygae/galgebra/tree/master/examples) for more examples in Python. With some small changes, they can be ported to Julia.

## Installation

The package can be installed with the Julia package manager. From the Julia REPL, type `]` to enter the **Pkg REPL mode** (you'll see a prompt like `(v1.1) pkg>`) and run:

```
dev https://github.com/pygae/GAlgebra.jl.git
```

The installation process will take a while, because it will install [PyCall.jl](https://github.com/JuliaPy/PyCall.jl) and its dependencies and it will also use `pip` to install [GAlgebra](https://github.com/pygae/galgebra) if you haven't. If you wish to use a version other than the latest released version of GAlgebra that will be installed by `deps/build.jl`, you can install that via pip before/after the installation of `GAlgebra.jl`, e.g. `pip install -e <your local path to GAlgebra>` or `pip install https://github.com/pygae/galgebra/archive/master.zip`.

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

## Python Dependency

GAlgebra.jl requires the Python [galgebra](https://github.com/pygae/galgebra) package (version ≥0.6.0, <0.7.0).

**How it works:**

1. On `Pkg.build("GAlgebra")`, the build script checks if galgebra is already installed in the Python environment that PyCall uses.
2. If galgebra is found and its version is in the supported range, it is used as-is.
3. If galgebra is not found or its version is outside the supported range, galgebra 0.6.0 is installed via pip.

**Using your own Python environment:**

GAlgebra.jl respects the Python that PyCall is configured to use. To use a specific Python (e.g. from a virtualenv):

```bash
export PYTHON=/path/to/your/venv/bin/python
julia -e 'using Pkg; Pkg.build("PyCall"); Pkg.build("GAlgebra")'
```

**Diagnostics:**

Set the `GALGEBRA_DEBUG` environment variable to see which galgebra version and Python executable are in use:

```bash
GALGEBRA_DEBUG=1 julia -e 'using GAlgebra'
```

This will print the galgebra version and Python path at module load time.
