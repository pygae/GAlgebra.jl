export G

"""

    G(p::Integer, q::Integer, r::Integer)

A convenient method to instantiate a Geometric Algebra with `p` positive, `q` negative and `r` zero dimensions.

`q` and `r` defaults to zero if omitted.

The basis will be named with prefix `e` and indexed starting from 1.

For example:

```julia
# Basic 
Hyper   = G(1)       # Hyperbolic numbers. 
ℂ       = G(0,1)     # Complex numbers.
Dual    = G(0,0,1)   # Dual numbers.
ℍ       = G(0,2)     # Quaternions.

# Clifford
Cl2 =       G(2)     # Clifford algebra for 2D vector space.
Cl3 =       G(3)     # Clifford algebra for 3D vector space.
Spacetime = G(1,3)   # Clifford algebra for timespace vectors.

# Geometric
PGA2D = G(2,0,1)     # Projective Euclidean 2D plane. (dual)
PGA3D = G(3,0,1)     # Projective Euclidean 3D space. (dual)
CGA2D = G(3,1)       # conformal 2D space. 
CGA3D = G(4,1)       # Conformal 3D space. 

# High-Dimensional GA
DCGA3D = G(6,2)      # Double Conformal 3D Space.
TCGA3D = G(9,3)      # Tripple Conformal 3D Space.
DCGSTA = G(4,8)      # Double Conformal Geometric Space Time Algebra.
QCGA   = G(9,6)      # Quadric Conformal Geometric Algebra.  
```

To instantiate a Geometric Algebra with more parameters, use `galgebra.ga.Ga` instead.

For example:

```julia
import SymPy: sympy
using GAlgebra

Ga = galgebra.ga.Ga

g3d = Ga("e*x|y|z")

(r, th, phi) = coords = sympy.symbols("r theta phi")
s3d = Ga("e_r e_theta e_phi", g=[1 0 0; 0 r^2 0; 0 0 r^2 * sympy.sin(th)^2], coords=coords, norm=true)
(er, eth, ephi) = s3d.mv()
```

Please also consult the documentation of [GAlgebra](https://github.com/pygae/galgebra).

"""
function G(p::Integer, q::Integer, r::Integer)
    total = p + q + r
    basis = "e*" * join([string(i) for i in 1:total], "|")
    metric = [fill(1, p) ; fill(-1, q) ; fill(0, r)]
    galgebra.ga.Ga(basis, g=metric)
end

G(p::Integer, q::Integer) = G(p, q, 0)

G(p::Integer) = G(p, 0, 0)

