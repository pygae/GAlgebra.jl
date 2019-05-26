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
    @test u ⋅ v == u | v == (u < v) == (u > v)
    @test v ⨼ B == (v < B)
    @test v ⨽ B == (v > B)

    @test v * v == (v * v).scalar()
    @test v * B == v ⋅ B + v ∧ B == v ⨼ B + v ∧ B
    @test u ⋅ v == u ⨼ v == u ⨽ v == u ⨰ v
    @test u ∧ v == u ⨱ v
    
    @test v == v.get_grade(1)
    @test G == G.get_grade(2)

    for r = 0:4
        @test (A + B).get_grade(r) == A.get_grade(r) + B.get_grade(r)
        @test (λ * A).get_grade(r) == (A * λ).get_grade(r) == λ * A.get_grade(r)

        Ar = A.get_grade(r)

        @test v ⨼ Ar == (v * Ar - (-1)^r * Ar * v) / 2
        @test Ar ⨽ v == (Ar * v - (-1)^r * v * Ar) / 2 == (-1)^(r-1) * (v ⨼ Ar)
        @test v ∧ Ar == (v * Ar + (-1)^r * Ar * v) / 2
        @test Ar ∧ v == (Ar * v + (-1)^r * v * Ar) / 2 == (-1)^r * (v ∧ Ar)

        @test v ⨼ Ar == (v * Ar).get_grade(r-1)
        @test v ∧ Ar == (v * Ar).get_grade(r+1)
        @test Ar ⨽ v == (Ar * v).get_grade(r-1)
        @test Ar ∧ v == (Ar * v).get_grade(r+1)
        @test v * Ar == v ⨼ Ar + v ∧ Ar
        @test Ar * v == Ar ⨽ v + Ar ∧ v

        Br = B.get_grade(r)
        Ar ⨼ Br == Ar ⨽ Br == (Ar * Br).scalar()

        for s = 0:4
            @test A.get_grade(r).get_grade(s) == (if r == s; A.get_grade(r) else 0 end)

            Bs = B.get_grade(s)
            ArBs = Ar * Bs
            @test ArBs == sum([ArBs.get_grade(abs(r - s) + 2 * j) for j=0:min(r, s)])

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

            for t in 0:4
                Ct = C.get_grade(t)

                Ar ∧ (Bs ∧ Ct) == (Ar * Bs * Ct).get_grade(r + s + t)
            end
        end
    end

    @test A == sum([A.get_grade(r) for r in 0:4])
    @test A.get_grade(-3) == 0

    @test v ⨼ A.even() == - (A.even() ⨽ v)
    @test v ⨼ A.odd() == A.odd() ⨽ v
    @test v ∧ A.even() == A.even() ∧ v
    @test v ∧ A.odd() == - (A.odd() ∧ v)

    @test A ⨼ B == sum([sum([(A.get_grade(r) * B.get_grade(s)).get_grade(s - r) for r in 0:4]) for s in 0:4])
    @test A ⨽ B == sum([sum([(A.get_grade(r) * B.get_grade(s)).get_grade(r - s) for r in 0:4]) for s in 0:4])
    @test A ∧ B == sum([sum([(A.get_grade(r) * B.get_grade(s)).get_grade(r + s) for r in 0:4]) for s in 0:4])

    @test (A ∧ B) ∧ C == A ∧ (B ∧ C) == A ∧ B ∧ C
    @test A ⨼ (B ⨽ C) == (A ⨼ B) ⨽ C
    @test A ⨼ (B ⨼ C) == (A ∧ B) ⨼ C
    @test A ⨽ (B ∧ C) == (A ⨽ B) ⨽ C
    @test (A ∧ B) ⨼ C == A ⨼ (B ⨼ C)

    @test u ∧ A ∧ v == - v ∧ A ∧ u

    @test A ⨰ B == A << B == (A * B + B * A) / 2
    @test A ⨱ B == A >> B == (A * B - B * A) / 2

    @test A ⨰ B == B ⨰ A
    @test A ⨱ B == - B ⨱ A

    @test A ⊛ B == B ⊛ A
    @test A ⊛ B == A.dual() ⊛ B.dual()
    @test A ⊛ B == A.rev() ⊛ B.rev()
    @test A ⊛ (B * C) == (B.rev() * A) ⊛ C
    @test A ⊛ (B ⨽ C) == (B.rev() ⨽ A) ⊛ C
    @test A ⊛ (B ⨼ C) == (B.rev() ∧ A) ⊛ C
    @test A ⊛ (B ∧ C) == (B.rev() ⨼ A) ⊛ C

    @test A * B ⋅ C ∧ D == ((A * B) ⋅ C) ∧ D

    @test u.dual() == u * V.I()
    @test v.project_in_blade(u) == (v ⋅ u) / u == (v ⨼ u) ⨼ u.inv()
    @test v.project_in_blade(w) + u.project_in_blade(w) == (u + v).project_in_blade(w)

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
