using Test
include("diagonal_sbp_D2.jl")

let
  bm = [1 4 6 8 11] # boundary stencil size
  N = 100
  for p = 2:2:10
    (D2, BS, HI, H, r) = diagonal_sbp_D2(p, N)

    # Check D2
    @test norm(sum(D2, dims=2)) < N * 1e-12
    @test norm(D2 * r) < N * 1e-12

    # Check with boundary
    for k = 2:p/2 + 1
      @test isapprox(D2 * (r.^k), k * (k-1) * r.^(k-2))
    end

    # Check without boundary
    b = bm[div(p, 2)]+1
    for k = 2:p+1
      x = D2 * (r.^k)
      y = k * (k-1) * r.^(k-2)
      @test isapprox(x[b:N+1-b], y[b:N+1-b])
    end

    # Check boundary stencil
    @test norm(sum(BS, dims=2)) < N * 1e-12
    for k = 1:div(p,2)+1
      x = BS * (r.^k)
      y = k * r.^(k-1)
      @test isapprox(x[1], -y[1])
      @test isapprox(x[end], y[end])
    end
  end
end
