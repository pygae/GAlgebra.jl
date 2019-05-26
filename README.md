# GAlgebra.jl

Julia interface to [GAlgebra](https://github.com/pygae/galgebra), a symbolic Geometric Algebra/Calculus package for SymPy.

## Development Status

Very early. But it already works and has many tests.

| **Build Status**                                                                                |
|:-----------------------------------------------------------------------------------------------:|
| [![][travis-img]][travis-url]                                   [![][codecov-img]][codecov-url] |

[travis-img]: https://travis-ci.com/pygae/GAlgebra.jl.svg?branch=master
[travis-url]: https://travis-ci.com/pygae/GAlgebra.jl
[codecov-img]: https://img.shields.io/codecov/c/github/pygae/GAlgebra.jl.svg
[codecov-url]: https://codecov.io/gh/pygae/GAlgebra.jl

## Getting Started

GAlgebra.jl itself doesn't depend on [SymPy.jl](https://github.com/JuliaPy/SymPy.jl), but is designed to work with it.

After installing SymPy.jl and GAlgebra.jl (see below for instructions), you may experiment with GAlgebra.jl just like in the Python version of GAlgebra (though there're some syntax differences like `True`/`true`, `'`/`"` etc. between Python and Julia).

For example, you may:

```julia
import SymPy: symbols, sympy, Sym
using GAlgebra

# In console, uncomment to enable colored printing with ANSI escape sequences 
# galgebra.printer.Eprint()
# In Jupyter, uncomment to enable LaTeX printing with MathJax
# galgebra.printer.Format()

(x, y, z) = xyz = symbols("x,y,z",real=true)
(o3d, ex, ey, ez) = galgebra.ga.Ga.build("e", g=[1, 1, 1], coords=xyz)

v = o3d.mv("v", "vector")
A = o3d.mv("A", "mv")

# Wedge product: \wedge
v ∧ A
# Hestene's inner product: \cdot
v ⋅ A
# Left contraction: \intprod
v ⨼ A
# Right contraction: \intprodr
v ⨽ A
# Scalar product: \circledast A ⊛ B = <A B†>
v ⊛ A
# Anti-comutator product: \dottimes  A⨰B = (AB+BA)/2
v ⨰ A
# Comutator product: \timesbar  A⨱B = (AB-BA)/2
v ⨱ A
```

Note: enter unicode symbols like `∧` with corresponding LaTeX like `\wedge` by [Tab completion](https://pkg.julialang.org/docs/julia/THl1k/1.1.0/manual/unicode-input.html).

So far only `galgebra.ga.Ga` and `galgebra.mv.Mv` have been verified to work in Julia, see [tests](https://github.com/pygae/GAlgebra.jl/tree/master/test/runtests.jl). The tests verified many identities in Linear Algebra and Geometric Algebra.

See [examples of GAlgebra](https://github.com/pygae/galgebra/tree/15-print-pow/examples) for more examples in Python. With some small changes, they can be ported to Julia.

## Installation

The package can be installed with the Julia package manager. From the Julia REPL, type `]` to enter the **Pkg REPL mode** (you'll see a prompt like `(v1.1) pkg>`) and run:

```
dev https://github.com/pygae/GAlgebra.jl.git
```

The installation process will take a while, because it will install [PyCall.jl](https://github.com/JuliaPy/PyCall.jl) and its dependencies and it will also use `pip` to install [GAlgebra](https://github.com/pygae/galgebra) (it currently requires `sympy==1.3` so this prerequisite is also ensured).

Now you may run GAlgebra.jl tests in Julia **Pkg REPL mode**: 

```
test GAlgebra
```

At the first time it will take a while, because SymPy.jl is specified as a test dependency so it'll be installed.

Then you'll see something like:

```
   Testing GAlgebra
 Resolving package versions...
Test Summary: | Pass  Total
GAlgebra.jl   |  300    300
   Testing GAlgebra tests passed
```

Hint: To get back to the Julia REPL please press backspace, see [Pkg doc](https://docs.julialang.org/en/v1/stdlib/Pkg/index.html) to learn more.