using LinearAlgebra

include("diagonal_sbp.jl")

let
  N = 20
  λ = [mod(i, 2) + 1 for i = 0:N]
  for p in (2, 4, 6)
    (_, _, _, _, _, A, _) = variable_diagonal_sbp_D2(p, N, λ; xc = (-1,1))
    (D1, _, H, _) = diagonal_sbp_D1(p, N; xc = (-1,1))
    R = A - D1' * H * Diagonal(λ) * D1
    eigenvalues = eigen(Matrix(R)).values
    println("SBP order = $p")
    @show extrema(real.(eigenvalues))
    @show extrema(imag.(eigenvalues))
    println()
  end
end
