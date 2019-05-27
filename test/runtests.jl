using PyCall
import SymPy: symbols, sympy, Sym
using GAlgebra
using Test

py"""
def vector(ga, components):
    bases = ga.mv()
    return sum([components[i] * e for i, e in enumerate(bases)])
"""
const vector = py"vector"

@testset "GAlgebra.jl" begin
    # Write your own tests here.
    (x, y, z) = xyz = symbols("x,y,z",real=true)
    (o3d, ex, ey, ez) = galgebra.ga.Ga.build("e", g=[1, 1, 1], coords=xyz)

    V = o3d
    dimV = range(0, stop=V.n)
    I = V.I()

    α = V.mv("α", "scalar")
    β = V.mv("β", "scalar")
    γ = V.mv("γ", "scalar")
    λ = V.mv("λ", "scalar")

    u = V.mv("u", "vector")
    v = V.mv("v", "vector")
    w = V.mv("w", "vector")

    A = V.mv("A", "mv")
    B = V.mv("B", "mv")
    C = V.mv("C", "mv")
    D = V.mv("D", "mv")

    G = V.mv("G", "grade", 2)
    R = V.mv("R", "spinor")

    # The following tests verified implementation correctness per definition

    @test u ⋅ v == u | v == (u < v) == (u > v) == u ⨼ v == u ⨽ v == u ⨰ v
    @test u ∧ v == u ⨱ v
    @test v ⨼ B == (v < B)
    @test v ⨽ B == (v > B)
    @test A ⨰ B == A << B == (A * B + B * A) / 2
    # @test A ×̄ B == A ⨰ B
    @test A ⨱ B == A >> B == (A * B - B * A) / 2
    @test A ⊛ B == A % B

    @test u × v == -I * (u ∧ v)
    @test_throws PyCall.PyError A × B

    @test -v == -1 * v
    @test abs(v) == v.norm()
    @test abs(R) == R.norm()
    @test ~A == A.rev()
    @test A' == adjoint(A) == A.dual() == A * I # Ga.dual_mode_value is default to "I+"
    @test v^-1 == inv(v) == v.inv() == v / (v^2)
    @test R^-1 == inv(R) == R.inv() == R / (R^2)

    @test v^-2 == (v^2).inv()
    @test R^-2 == (R^2).inv()
    @test v^0 == 1
    @test v^2 == v*v

    @test (v)⁻¹ == v^-1
    @test (R)⁻¹ == R^-1
    @test (A)⁻ == involution(A) == A.even() - A.odd()
    @test (A)ᵀ == A.rev()
    @test (A)ǂ == involution(A).rev()
    @test proj(u, v) == v.project_in_blade(u)
    @test refl(u, v) == v.reflect_in_blade(u)
    @test rot(u ∧ v, A) == A.rotate_multivector(u ∧ v)
    @test exp(u ∧ v) == (u ∧ v).exp()

    for r ∈ dimV
        A[r] == A.grade(r) == A.get_grade(r)
    end

    # The following tests verified many identities in Linear Algebra

    @test v + w == w + v
    @test (u + v) + w == u + (v + w)
    @test v + 0 == v
    @test 0 * v == 0
    @test 1 * v == v
    @test α * (β * v) == (α * β) * v
    @test α * (v + w) == α * v + α * w
    @test (α + β) * v == α * v + β * v
    @test v + (-1) * v == 0

    𝑶 = vector(V, [0, 0, 0])
    @test α * 𝑶 == 𝑶
    @test (-α) * v == α * (-v) == -α * v

    # The following tests verified many identities in https://arxiv.org/abs/1205.5935

    @test v * v == (v * v).scalar()
    @test v * B == v ⋅ B + v ∧ B == v ⨼ B + v ∧ B

    @test u ∧ (v + λ * u) == u ∧ v
    
    @test v == v[1]
    @test G == G[2]

    for r ∈ dimV
        @test (A + B)[r] == A[r] + B[r]
        @test (λ * A)[r] == (A * λ)[r] == λ * A[r]

        Ar = A[r]

        @test v ⨼ Ar == (v * Ar - (-1)^r * Ar * v) / 2
        @test Ar ⨽ v == (Ar * v - (-1)^r * v * Ar) / 2 == (-1)^(r-1) * (v ⨼ Ar)
        @test v ∧ Ar == (v * Ar + (-1)^r * Ar * v) / 2
        @test Ar ∧ v == (Ar * v + (-1)^r * v * Ar) / 2 == (-1)^r * (v ∧ Ar)

        @test v ⨼ Ar == (v * Ar)[r-1]
        @test v ∧ Ar == (v * Ar)[r+1]
        @test Ar ⨽ v == (Ar * v)[r-1]
        @test Ar ∧ v == (Ar * v)[r+1]
        @test v * Ar == v ⨼ Ar + v ∧ Ar
        @test Ar * v == Ar ⨽ v + Ar ∧ v

        Br = B[r]
        Ar ⨼ Br == Ar ⨽ Br == (Ar * Br).scalar()

        for s ∈ dimV
            @test A[r][s] == (if r == s; A[r] else 0 end)

            Bs = B[s]
            ArBs = Ar * Bs

            @test ArBs == sum([ArBs[abs(r - s) + 2j] for j=0:min(r, s)])                      # A.4.1
            @test Ar ⨼ Bs == (-1)^(r * (s - 1)) * Bs ⨽ Ar                                    # A.4.10
            @test Ar ∧ Bs == (-1)^(r * s) * Bs ∧ Ar                                          # A.4.11

            for j ∈ dimV
                @test ArBs[r + s - 2j] == (-1)^(r * s - j) * (B[s] * A[r])[r + s - 2j]        # A.4.2
            end

            @test v ⨼ ArBs == (v * ArBs - (-1)^(r+s) * ArBs * v)/2 ==
                (v ⨼ Ar) * Bs + (-1)^r * Ar * (v ⨼ Bs) == 
                (v ∧ Ar) * Bs - (-1)^r * Ar * (v ∧ Bs)
            @test v ∧ ArBs == (v * ArBs + (-1)^(r+s) * ArBs * v)/2 ==
                (v ∧ Ar) * Bs - (-1)^r * Ar * (v ⨼ Bs) ==
                (v ⨼ Ar) * Bs + (-1)^r * Ar * (v ∧ Bs)

            @test v ⨼ (Ar ∧ Bs) == (v ⨼ Ar) ∧ Bs + (-1)^r * Ar ∧ (v ⨼ Bs)
            @test v ∧ (Ar ⨽ Bs) == (v ∧ Ar) ⨽ Bs - (-1)^r * Ar ⨽ (v ⨼ Bs)
            @test v ∧ (Ar ⨼ Bs) == (v ⨼ Ar) ⨼ Bs + (-1)^r * Ar ⨼ (v ∧ Bs)

            if r > s
                @test Ar ⨼ Bs == Bs ⨽ Ar == 0
            end

            for t ∈ dimV
                Ct = C[t]

                Ar ∧ (Bs ∧ Ct) == (Ar * Bs * Ct)[r + s + t]
            end
        end
    end

    @test A == sum([A[r] for r ∈ dimV])
    @test A[-3] == 0

    @test v ⨼ A == (v * A - (A)⁻ * v)/2                                                      # A.4.13
    @test v ∧ A == (v * A + (A)⁻ * v)/2                                                      # A.4.14
    @test A ⨽ v == - v ⨼ (A)⁻                                                                # A.4.15
    @test A ∧ v == v ∧ (A)⁻                                                                  # A.4.16

    @test v ⨼ (A * B) == (v ⨼ A) * B + (A)⁻ * (v ⨼ B) == (v ∧ A) * B - (A)⁻ * (v ∧ B)       # A.4.18-19
    @test v ∧ (A * B) == (v ∧ A) * B - (A)⁻ * (v ⨼ B) == (v ⨼ A) * B + (A)⁻ * (v ∧ B)       # A.4.20-21
    @test v ⨼ (A ∧ B) == (v ⨼ A) ∧ B + (A)⁻ ∧ (v ⨼ B)                                       # A.4.22
    @test v ∧ (A ⨽ B) == (v ∧ A) ⨽ B - (A)⁻ ⨽ (v ⨼ B)                                       # A.4.23
    @test v ∧ (A ⨼ B) == (v ⨼ A) ⨼ B + (A)⁻ ⨼ (v ∧ B)                                       # A.4.24

    @test v ⨼ A.even() == - (A.even() ⨽ v)
    @test v ⨼ A.odd() == A.odd() ⨽ v
    @test v ∧ A.even() == A.even() ∧ v
    @test v ∧ A.odd() == - (A.odd() ∧ v)

    @test (A * B).scalar() == (B * A).scalar() == (~A * ~B).scalar() == 
        ((A)⁻ * (B)⁻).scalar() == ((A)ᵀ * (B)ᵀ).scalar() == ((A)ǂ * (B)ǂ).scalar()            # A.4.3-6

    @test A ⨼ B == sum([sum([(A[r] * B[s])[s - r] for r ∈ dimV]) for s ∈ dimV])             # A.4.7
    @test A ⨽ B == sum([sum([(A[r] * B[s])[r - s] for r ∈ dimV]) for s ∈ dimV])             # A.4.8
    @test A ∧ B == sum([sum([(A[r] * B[s])[r + s] for r ∈ dimV]) for s ∈ dimV])             # A.4.9

    @test (A ∧ B) ∧ C == A ∧ (B ∧ C) == A ∧ B ∧ C                                           # A.4.28
    @test A ⨼ (B ⨽ C) == (A ⨼ B) ⨽ C                                                        # A.4.29
    @test A ⨼ (B ⨼ C) == (A ∧ B) ⨼ C                                                        # A.4.30
    @test A ⨽ (B ∧ C) == (A ⨽ B) ⨽ C                                                        # A.4.31
    @test (A ∧ B) ⨼ C == A ⨼ (B ⨼ C)

    @test u ∧ A ∧ v == - v ∧ A ∧ u                                                           # A.4.17

    @test A * B == A ⨱ B + A ⨰ B
    @test A ⨰ B == B ⨰ A
    @test A ⨱ B == - B ⨱ A

    @test A ⊛ B == B ⊛ A
    @test A ⊛ B == A' ⊛ B' == A.dual() ⊛ B.dual()
    @test A ⊛ B == ~A ⊛ ~B == A.rev() ⊛ B.rev()
    @test A ⊛ (B * C) == (~B * A) ⊛ C
    @test A ⊛ (B ⨽ C) == (~B ⨽ A) ⊛ C
    @test A ⊛ (B ⨼ C) == (~B ∧ A) ⊛ C
    @test A ⊛ (B ∧ C) == (~B ⨼ A) ⊛ C

    @test A * B ⋅ C ∧ D == ((A * B) ⋅ C) ∧ D

    @test u.dual() == u * V.I()
    @test proj(u, v) == (v ⋅ u) / u == (v ⨼ u) ⨼ u.inv()
    @test proj(w, v) + proj(w, u) == proj(w, u + v)

    Vr = u ∧ v
    @test proj(Vr, B) == B ⨼ Vr * (Vr)⁻¹ == (B ⨼ Vr) ⨼ (Vr)⁻¹                               # A.4.34
    # TODO this is failing for now
    # @test refl(Vr, B) == B ∧ Vr * (Vr)⁻¹ == (B ∧ Vr) ⨽ (Vr)⁻¹                               # A.4.35

    # The following tests verified interoperability with numeric and symbolic numbers

    uu = vector(V, [1, 2, 3])
    vv = vector(V, [4, 5, 6])
    ww = vector(V, [5, 6, 7])

    @test uu + vv == 5 * ex + 7 * ey + 9 * ez
    @test 7 * uu + 2 * ww == 17 * ex + 26 * ey + 35 * ez
    @test 7 * uu - 2 * ww == -3 * ex + 2 * ey + 7 * ez
    @test 3 * uu + 2 * vv + ww == 16 * ex + 22 * ey + 28 * ez
    @test (sympy.sqrt(2) * u + sympy.Rational(2, 3) * v) ⋅ ey == 
        sympy.sqrt(2) * (u ⋅ ey) + sympy.Rational(2, 3) * (v ⋅ ey)
end
