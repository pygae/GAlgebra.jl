# GAlgebra.jl API

```@index
```

```@autodocs
Modules = [GAlgebra]
Pages = ["ga.jl"]
```

## class Mv

```@autodocs
Modules = [GAlgebra]
Pages = ["mv.jl"]
```

```@docs
Base.:+(::Mv,::Mv)
Base.:-(::Mv,::Mv)
Base.:*(::Mv,::Mv)
Base.:/(::Mv,::Mv)
Base.:(==)(::Mv,::Mv)
GAlgebra.:≠(::Mv,::Mv)

GAlgebra.:∧(::Mv,::Mv)
GAlgebra.:⋅(::Mv,::Mv)
GAlgebra.:⨼(::Mv,::Mv)
GAlgebra.:⨽(::Mv,::Mv)

GAlgebra.:⊠(::Mv,::Mv)
GAlgebra.:⊙(::Mv,::Mv)
GAlgebra.:×(::Mv,::Mv)
GAlgebra.:⊛(::Mv,::Mv)

Base.:-(::Mv)

GAlgebra.norm(::Mv)
Base.inv(::Mv)
Base.:~(::Mv)
Base.adjoint(::Mv)
GAlgebra.involute(::Mv)
Base.conj(::Mv)

GAlgebra.proj(::Mv,::Mv)
GAlgebra.refl(::Mv,::Mv)
GAlgebra.rot(::Mv,::Mv)
Base.exp(::Mv)

GAlgebra.scalar(::Mv)
GAlgebra.even(::Mv)
GAlgebra.odd(::Mv)
```

