include("diagonal_sbp.jl")
using SparseArrays

⊗ = (A,B) -> kron(A, B)

function locoperator(p, Nx, Ny, corners::NTuple{4, Tuple{T, T}}) where T <: Number
  @assert (corners[1][1], corners[2][1]) == (corners[3][1], corners[4][1])
  @assert (corners[1][2], corners[3][2]) == (corners[2][2], corners[4][2])

  (D2x, BSx, HIx, Hx, rx) = diagonal_sbp_D2(p, Nx; xc = (corners[1][1],
                                                         corners[4][1]))
  (D2y, BSy, HIy, Hy, ry) = diagonal_sbp_D2(p, Ny; xc = (corners[1][2],
                                                         corners[4][2]))

  Ax= BSx- Hx* D2x
  Ay= BSy- Hy* D2y

  E0x = sparse([     1], [     1], [1], Nx + 1, Nx + 1)
  ENx = sparse([Nx + 1], [Nx + 1], [1], Nx + 1, Nx + 1)
  E0y = sparse([     1], [     1], [1], Ny + 1, Ny + 1)
  ENy = sparse([Ny + 1], [Ny + 1], [1], Ny + 1, Ny + 1)

  e0x = sparse([     1], [1], [1], Nx + 1, 1)
  eNx = sparse([Nx + 1], [1], [1], Nx + 1, 1)
  e0y = sparse([     1], [1], [1], Ny + 1, 1)
  eNy = sparse([Ny + 1], [1], [1], Ny + 1, 1)


  hx = rx[2] - rx[1]
  hy = ry[2] - ry[1]

  taux = 4/hx
  tauy = 4/hy
  M = (Hy ⊗ (Ax - BSx - BSx' + taux * (E0x + ENx))) +
      ((Ay - BSy - BSy' + tauy * (E0y + ENy)) ⊗ Hx)
  B0x = (Hy ⊗ (-BSx' * e0x + taux * e0x))
  BNx = (Hy ⊗ (-BSx' * eNx + taux * eNx))
  B0y = ((-BSy' * e0y + tauy * e0y) ⊗ Hx)
  BNy = ((-BSy' * eNy + tauy * eNy) ⊗ Hx)

  r = ones(Ny + 1) ⊗ rx
  s = ry ⊗ ones(Nx + 1)
  ((M, B0x, BNx, B0y, BNy), r, s, rx, ry, Hy ⊗ Hx)
end

function localsolve(M, g0x, gNx, g0y, gNy)
  localsolve(M..., g0x, gNx, g0y, gNy)
end

function localsolve(M, B0x, BNx, B0y, BNy, g0x, gNx, g0y, gNy)
  b = B0x * g0x + BNx * gNx + B0y * g0y + BNy * gNy
  M \ b
end

let
  # min size required for order p/2
  N0 = (2, 8, 12, 16, 21)

  # What order to run
  p = 4

  # number of levels run is based on size of this array
  ϵ = zeros(5)
  for lvl = 1:length(ϵ)
    # Setup the local operators
    Nx = 2^(lvl-1) * (N0[div(p,2)] + 0)
    Ny = 2^(lvl-1) * (N0[div(p,2)] + 2)
    (lop, x, y, x1d, y1d, H) =
      locoperator(p, Nx, Ny, ((0, 0), (1, 0), (0, 1), (1, 1)))

    # Make sure we're positive definite and symmetric
    if lvl == 2
      (λ, V) = eigen(Matrix(lop[1]))
      println((minimum(real.(λ)), maximum(real.(λ)), minimum(imag.(λ)),
               maximum(imag.(λ))))
      println(lop[1] ≈ lop[1]')
    end

    # Exact Solution
    uexact = (x,y) -> sin.(π * x) .* sinh.(π * y)

    # Dirichlet boundary conditions
    g0x = uexact(0, y1d)
    gNx = uexact(1, y1d)
    g0y = uexact(x1d, 0)
    gNy = uexact(x1d, 1)

    # Solve local problem
    u = localsolve(lop, g0x, gNx, g0y, gNy)

    # Check the error
    Δu = u - uexact(x, y)
    ϵ[lvl] = sqrt(Δu' * H * Δu)
  end
  # Check the rate
  println(ϵ)
  println((log.(ϵ[1:end-1]) - log.(ϵ[2:end])) / log(2))

  nothing
end
