using PyCall
import SymPy: symbols, sympy, Sym
using GAlgebra
using Test

py"""
def vector(ga, components):
    basis = ga.mv()
    return sum([components[i] * e for i, e in enumerate(basis)])
"""
const vector = py"vector"

py"""
def signature(ga):
    basis = ga.mv()
    signs = [e * e for e in basis]
    p, q, r = 0, 0, 0
    for sign in signs:
        p += 1 if sign == 1 else 0
        q += 1 if sign == -1 else 0
        r += 1 if sign == 0 else 0

    return (p, q, r)
"""
const signature = py"signature"

macro ctest(assertion::Expr, sep::QuoteNode, condition::Expr)
    if sep == :(:if)
        esc(quote
            if $condition
                @test $assertion
            else
                @test_broken $assertion
            end
        end)
    elseif sep == :(:unless)
        esc(quote
            if $condition
                @test_broken $assertion
            else
                @test $assertion
            end
        end)
    else
        throw(DomainError(sep, "sep can only be one of :if, :unless"))
    end    
end

@testset "GAlgebra.jl" begin
    # Basic 
    Hyper   = G(1)       # Hyperbolic numbers. 
    ℂ       = G(0,1)     # Complex numbers.
    Dual    = G(0,0,1)   # Dual numbers.
    ℍ       = G(0,2)     # Quaternions.

    # Clifford
    Cl2       = G(2)     # Clifford algebra for 2D vector space.
    Cl3       = G(3)     # Clifford algebra for 3D vector space.
    Spacetime = G(1,3)   # Clifford algebra for timespace vectors.

    # Geometric
    PGA2D = G(2,0,1)     # Projective Euclidean 2D plane. (dual)
    PGA3D = G(3,0,1)     # Projective Euclidean 3D space. (dual)
    CGA2D = G(3,1)       # Conformal 2D space. 
    CGA3D = G(4,1)       # Conformal 3D space. 

    # for V ∈ [PGA2D, PGA3D, CGA2D, CGA3D] # pass but take extremely long
    # Even PGA2D takes 189.979728 seconds (40.41 k allocations: 1.395 MiB)

    # for V ∈ [Hyper, Dual, ℂ] # for quick test any of p, q, r is 1

    for V ∈ [Cl2, Cl3, ℂ, ℍ, Hyper, Dual, Spacetime]
        sigV = signature(V)
        @time @testset "G$sigV" begin
            # @test_broken 1==1 # for triggering error
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

            R = V.mv("R", "spinor")

            # Precalculte AB and BA
            AB = A * B
            BA = B * A

            # The following tests verified implementation correctness per definition

            @test u ⋅ v == u | v == (u < v) == (u > v) == u ⨼ v == u ⨽ v == u ⨰ v
            @test u ∧ v == u ⨱ v
            @test v ⨼ B == (v < B)
            @test v ⨽ B == (v > B)
            @test A ⨰ B == A << B == (AB + BA) / 2
            # @test A ×̄ B == A ⨰ B
            @test A ⨱ B == A >> B == (AB - BA) / 2
            @test A ⊛ B == A % B

            @test abs(v) == norm(v) == v.norm()
            if V ∉ [Spacetime, PGA2D, PGA3D, CGA2D, CGA3D]
                @test abs(R) == norm(R) == R.norm()
            end

            @test ~A == A[:~] == rev(A) == A.rev()

            if V ∉ [Dual, PGA2D, PGA3D, CGA2D, CGA3D]
                @test A' == dual(A) == A.dual() == adjoint(A) == A * I # Ga.dual_mode_value is default to "I+"
                @test (v)⁻¹ == v[:⁻¹] == v^-1 == inv(v) == v.inv()
                @test v^-2 == (v^2).inv()
            end
            
            @test (A)ˣ == A[:*] == involute(A) == (A)₊ - (A)₋ == A[:+] - A[:-] == A.even() - A.odd()
            @test (A)ǂ == A[:ǂ] == conj(A) == involute(A).rev()   

            if V ∉ [Spacetime, PGA2D, PGA3D, CGA2D, CGA3D]
                @test R^-2 == (R^2).inv()
                @test (R)⁻¹ == R[:⁻¹] == R^-1 == inv(R) == R.inv()
                @test ((R)⁻¹)ˣ == ((R)ˣ)⁻¹
                @test ((R)⁻¹)ǂ == ((R)ǂ)⁻¹
            end

            if V ∈ [Cl2, Cl3]
                @test (v)⁻¹ == (~v) / norm(v)^2 == v / v^2 
                @test (R)⁻¹ == (~R) / norm(R)^2 == R / R^2
            end

            @test v^0 == 1
            @test v^2 == v*v
            
            @test ((A)ˣ)ˣ == ~(~A) == A[:~][:~] == ((A)ǂ)ǂ == A
            @test ~((A)ˣ) == (~A)ˣ

            if V ∉ [Dual]
                @test proj(u, v) == v.project_in_blade(u)
                @test refl(u, v) == v.reflect_in_blade(u)
            end
            
            @test rot(u ∧ v, A) == A.rotate_multivector(u ∧ v)
            @test exp(u ∧ v) == (u ∧ v).exp()

            @test typeof(scalar(A)) == Sym
            @test typeof(A[0]) == Mv
            @test scalar(A) == A.scalar() == A[0].obj
            @test (A)₊ == A[:+] == even(A) == A.even()
            @test (A)₋ == A[:-] == odd(A) == A.odd()

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
            @test -v == -1 * v

            𝑶 = vector(V, fill(0, V.n))
            @test α * 𝑶 == 𝑶
            @test (-α) * v == α * (-v) == -α * v

            # The following tests verified many identities in https://arxiv.org/abs/1205.5935

            @test v * v == (v * v).scalar()
            @test v * B == v ⋅ B + v ∧ B == v ⨼ B + v ∧ B

            @test u ∧ (v + λ * u) == u ∧ v
            
            @test v == v[1]
            if V.n >= 2
                G2 = V.mv("G2", "grade", 2)
                @test G2 == G2[2]
            end

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

            @test v ⨼ A == (v * A - (A)ˣ * v)/2                                                      # A.4.13
            @test v ∧ A == (v * A + (A)ˣ * v)/2                                                      # A.4.14
            @test A ⨽ v == - v ⨼ (A)ˣ                                                                # A.4.15
            @test A ∧ v == v ∧ (A)ˣ                                                                  # A.4.16

            @test v ⨼ (AB) == (v ⨼ A) * B + (A)ˣ * (v ⨼ B) == (v ∧ A) * B - (A)ˣ * (v ∧ B)       # A.4.18-19
            @test v ∧ (AB) == (v ∧ A) * B - (A)ˣ * (v ⨼ B) == (v ⨼ A) * B + (A)ˣ * (v ∧ B)       # A.4.20-21
            @test v ⨼ (A ∧ B) == (v ⨼ A) ∧ B + (A)ˣ ∧ (v ⨼ B)                                       # A.4.22
            @test v ∧ (A ⨽ B) == (v ∧ A) ⨽ B - (A)ˣ ⨽ (v ⨼ B)                                       # A.4.23
            @test v ∧ (A ⨼ B) == (v ⨼ A) ⨼ B + (A)ˣ ⨼ (v ∧ B)                                       # A.4.24

            @test v ⨼ A[:+] == - (A[:+] ⨽ v)
            @test v ⨼ A[:-] == A[:-] ⨽ v
            @test v ∧ A[:+] == A[:+] ∧ v
            @test v ∧ A[:-] == - (A[:-] ∧ v)

            @test (AB).scalar() == (BA).scalar() == (~A * ~B).scalar() == 
                ((A)ˣ * (B)ˣ).scalar() == ((A)ǂ * (B)ǂ).scalar()                                      # A.4.3-6

            @test A ⨼ B == sum([sum([(A[r] * B[s])[s - r] for r ∈ dimV]) for s ∈ dimV])             # A.4.7
            @test A ⨽ B == sum([sum([(A[r] * B[s])[r - s] for r ∈ dimV]) for s ∈ dimV])             # A.4.8
            @test A ∧ B == sum([sum([(A[r] * B[s])[r + s] for r ∈ dimV]) for s ∈ dimV])             # A.4.9

            @test (A ∧ B) ∧ C == A ∧ (B ∧ C) == A ∧ B ∧ C                                           # A.4.28
            @test A ⨼ (B ⨽ C) == (A ⨼ B) ⨽ C                                                        # A.4.29
            @test A ⨼ (B ⨼ C) == (A ∧ B) ⨼ C                                                        # A.4.30
            @test A ⨽ (B ∧ C) == (A ⨽ B) ⨽ C                                                        # A.4.31
            @test (A ∧ B) ⨼ C == A ⨼ (B ⨼ C)

            @test u ∧ A ∧ v == - v ∧ A ∧ u                                                           # A.4.17

            @test AB == A ⨱ B + A ⨰ B
            @test A ⨰ B == B ⨰ A
            @test A ⨱ B == - B ⨱ A

            @test A ⊛ B == B ⊛ A
            if V ∉ [Spacetime, ℂ, ℍ, Dual, PGA2D, PGA3D, CGA2D, CGA3D]
                @test A ⊛ B == A' ⊛ B' == A.dual() ⊛ B.dual()
            end      
            @test A ⊛ B == ~A ⊛ ~B == A.rev() ⊛ B.rev()
            @test A ⊛ (B * C) == (~B * A) ⊛ C
            @test A ⊛ (B ⨽ C) == (~B ⨽ A) ⊛ C
            @test A ⊛ (B ⨼ C) == (~B ∧ A) ⊛ C
            @test A ⊛ (B ∧ C) == (~B ⨼ A) ⊛ C

            @test AB ⋅ C ∧ D == ((AB) ⋅ C) ∧ D

            if V ∉ [Dual, PGA2D, PGA3D, CGA2D, CGA3D]
                @test u.dual() == u * V.I()
                @test proj(u, v) == (v ⋅ u) / u == (v ⨼ u) ⨼ u.inv()
                @test proj(w, v) + proj(w, u) == proj(w, u + v)
            end

            if V == Cl3
                @test u × v == -I * (u ∧ v)
                @test_throws PyCall.PyError A × B

                Vr = u ∧ v
                @test proj(Vr, B) == B ⨼ Vr * (Vr)⁻¹ == (B ⨼ Vr) ⨼ (Vr)⁻¹                               # A.4.34
                # TODO this is failing for now
                @test_broken refl(Vr, B) == B ∧ Vr * (Vr)⁻¹ == (B ∧ Vr) ⨽ (Vr)⁻¹                               # A.4.35

                # The following tests verified interoperability with numeric and symbolic numbers
                (ex, ey, ez) = V.mv()

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
        end

        # The following takes forever

        # High-Dimensional GA
        # DCGA3D = G(6,2)      # Double Conformal 3D Space.
        # TCGA3D = G(9,3)      # Tripple Conformal 3D Space.
        # DCGSTA = G(4,8)      # Double Conformal Geometric Space Time Algebra.
        # QCGA   = G(9,6)      # Quadric Conformal Geometric Algebra. 
    end
end