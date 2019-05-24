# GAlgebra.jl

Julia interface to [GAlgebra](https://github.com/pygae/galgebra), a symbolic Geometric Algebra/Calculus package for SymPy.

## Development Status

Very early. But it already works and has some unit tests.

## Getting Started

GAlgebra.jl itself doesn't depend on [SymPy.jl](https://github.com/JuliaPy/SymPy.jl), but is designed to work with it.

After installing SymPy.jl and GAlgebra.jl (see below for instructions), you may experiment with GAlgebra.jl just like in Python (though there're some syntax differences like `True`/`true`, `'`/`"` etc. between Python and Julia).

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
v | A
# Left contraction: \intprod
v ⨼ A
v < A
# Right contraction: \intprodr
v ⨽ A
v > A
# Anti-comutator product: A<<B = (AB+BA)/2
v << A
# Comutator product: \timesbar  A⨱B = (AB-BA)/2
v ⨱ A
v >> A
```

Note: enter unicode symbols like `∧` with corresponding LaTeX like `\wedge` by [Tab completion](https://pkg.julialang.org/docs/julia/THl1k/1.1.0/manual/unicode-input.html).

So far only `galgebra.ga.Ga` and `galgebra.mv.Mv` have been verified to work in Julia, see [tests](https://github.com/pygae/GAlgebra.jl/tree/master/test/runtests.jl).

See [examples of GAlgebra](https://github.com/pygae/galgebra/tree/15-print-pow/examples) for more examples in Python. With some small changes, they can be ported to Julia.

## Installation

Ideally, you should be able to install GAlgebra.jl by the following command in Julia **Pkg** REPL:

```
add GAlgebra
```

But it doesn't work like this yet.

For the time being, you need to do the following 5 steps:

Step 1: install PyCall.jl in Julia **Pkg** REPL:

```
add PyCall
```

Hint: pressing `]` will enter Pkg REPL mode (you'll see a prompt like `(v1.1) pkg>`), to get back to the Julia REPL please press backspace, see [Pkg doc](https://docs.julialang.org/en/v1/stdlib/Pkg/index.html) to learn more.

Step 2：install from the latest source code of GAlgebra for the Python environment used by PyCall.jl in **Julia** REPL:

```julia
using PyCall
run(PyCall.python_cmd(`-m pip install sympy==1.3 -e git+https://github.com/pygae/galgebra.git#egg=galgebra`))
```

Step 3：install SymPy.jl in Julia **Pkg** REPL:

```
add SymPy
```

Step 4: install GAlgebra.jl as a developing package in Julia **Pkg** REPL: 

```
dev https://github.com/pygae/GAlgebra.jl.git
```

Step 5: run GAlgebra.jl unit tests in Julia **Pkg** REPL: 

```
test GAlgebra
```

Then you'll see something like:

```
   Testing GAlgebra
 Resolving package versions...
Test Summary: | Pass  Total
GAlgebra.jl   |   36     36
   Testing GAlgebra tests passed
```
