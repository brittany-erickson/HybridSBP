do_plotting = true
if VERSION < v"0.6.999999"
  macro isdefined(s::Symbol)
    return isdefined(s)
  end
end
if !@isdefined mtime_global_curved
  mtime_global_curved = 0
end
if mtime_global_curved < mtime("global_curved.jl")
  println("including global_curved")
  include("global_curved.jl")
  mtime_global_curved = mtime("global_curved.jl")
end
if VERSION < v"0.6.999999"
  macro isdefined(s::Symbol)
    return isdefined(s)
  end
end
if !@isdefined mtime_global_curved
  mtime_global_curved = 0
end
if VERSION <= v"0.6.999999"
  ldiv! = A_ldiv_B!
  cholesky = cholfact
  using CholmodSolve2
end

let
  TOP_MAIN = 7
  BOTTOM_MAIN = 8
  BRANCH = 9

  SBPp   = 4 # SBP interior order

  (verts, EToV, EToF, FToB, EToBlock) = read_inp_2d("meshes/branch.inp")
  Lx = Ly = 10
  H = 0.4 * Lx
  verts = Lx*verts
  N1 = N0 = 13


  # number of elements and faces
  (nelems, nfaces) = (size(EToV, 2), size(FToB, 1))
  @show (nelems, nfaces)

  EToN0 = zeros(Int64, 2, nelems)
  EToN0[1, :] .= N0
  EToN0[2, :] .= N1

  @assert typeof(EToV) == Array{Int, 2} && size(EToV) == (4, nelems)
  @assert typeof(EToF) == Array{Int, 2} && size(EToF) == (4, nelems)
  @assert maximum(maximum(EToF)) == nfaces

  # Determine secondary arrays
  # FToE : Unique Global Face to Element Number
  # FToLF: Unique Global Face to Element local face number
  # EToO : Element to Unique Global Faces Orientation
  # EToS : Element to Unique Global Face Side
  (FToE, FToLF, EToO, EToS) = connectivityarrays(EToV, EToF)

  # Exact solution
  u1ex(x,y) = cos.(x) .* cos.(π * y / Ly) .- (H .* x - (x.^2)/2)
  u1ex_x(x,y) = -sin.(x) .* cos.(π * y / Ly) .- (H .- x)
  u1ex_y(x,y) = -(π / Ly) * cos.(x) .* sin.(π * y / Ly)
  u1ex_xx(x,y) = -cos.(x) .* cos.(π * y / Ly) .+ 1
  u1ex_yy(x,y) = -(π/Ly)^2 .* cos.(x) .* cos.(π * y / Ly)

  u2ex(x,y) = 1 .+ sin.(x.^2 .+ (y .+ H).^2) .+ (y.^2)/2 .- (H * x .+ (x.^2)/2)
  u2ex_x(x,y) = -H .- x .+ 2 .* x .* cos.(x.^2 .+ (H .+ y).^2)
  u2ex_y(x,y) = y .+ 2 .* (H .+ y) .* cos.(x.^2 .+ (H .+ y).^2)
  u2ex_xx(x,y) = 2 * (cos.(x.^2 .+ (y .+ H).^2) .- 2 * x.^2 .* sin.(x.^2 + (y .+ H).^2)) .- 1
  u2ex_yy(x,y) = 2 * (cos.(x.^2 .+ (H .+ y).^2) .- 2 .* (H .+ y).^2 .* sin.(x.^2 + (H .+ y).^2)) .+ 1

  u3ex(x,y) = cos.(x.^2 + (y .+ H).^2) + (y.^2)/2 .- (H * x + (1/2) * x.^2)
  u3ex_x(x,y) = -H .- x .- 2 .* x .* sin.(x.^2 .+ (H .+ y).^2)
  u3ex_y(x,y) = y .- 2 .* (H .+ y) .* sin.(x.^2 .+ (H .+ y).^2)
  u3ex_xx(x,y) = -1 .- 4 .* x.^2 .* cos.(x.^2 .+ (H .+ y).^2) .- 2 .* sin.(x.^2 .+ (H .+ y).^2)
  u3ex_yy(x,y) = 1 .- 4 .* (H .+ y).^2 .* cos.(x.^2 .+ (H .+ y).^2) .- 2 .* sin.(x.^2 .+ (H .+ y).^2)

  vex(x,y,e) = begin
    if EToBlock[e] == 1
      return u1ex(x,y)
    elseif EToBlock[e] == 2
      return u2ex(x,y)
    elseif EToBlock[e] == 3
      return u3ex(x,y)
    else
      error("invalid block")
    end
  end
  vex_x(x,y,e) = begin
    if EToBlock[e] == 1
      return u1ex_x(x,y)
    elseif EToBlock[e] == 2
      return u2ex_x(x,y)
    elseif EToBlock[e] == 3
      return u3ex_x(x,y)
    else
      error("invalid block")
    end
  end
  vex_y(x,y,e) = begin
    if EToBlock[e] == 1
      return u1ex_y(x,y)
    elseif EToBlock[e] == 2
      return u2ex_y(x,y)
    elseif EToBlock[e] == 3
      return u3ex_y(x,y)
    else
      error("invalid block")
    end
  end
  fex(x,y,e) = begin
    if EToBlock[e] == 1
      return u1ex_xx(x,y) + u1ex_yy(x,y)
    elseif EToBlock[e] == 2
      return u2ex_xx(x,y) + u2ex_yy(x,y)
    elseif EToBlock[e] == 3
      return u3ex_xx(x,y) + u3ex_yy(x,y)
    else
      error("invalid block")
    end
  end

  ϵ = zeros(4)
  for lvl = 1:length(ϵ)
    # Set up the local grid dimensions
    Nr = EToN0[1, :] * (2^(lvl-1))
    Ns = EToN0[2, :] * (2^(lvl-1))

    #{{{ Build the local volume operators
    # Dictionary to store the operators
    OPTYPE = typeof(locoperator(2, 8, 8, (r,s)->r, (r,s)->s))
    lop = Dict{Int64, OPTYPE}()
    for e = 1:nelems
      # @show (e, nelems)
      # Get the element corners
      (x1, x2, x3, x4) = verts[1, EToV[:, e]]
      (y1, y2, y3, y4) = verts[2, EToV[:, e]]

      xt = (r,s)->transfinite_blend(x1, x2, x3, x4, r, s)
      yt = (r,s)->transfinite_blend(y1, y2, y3, y4, r, s)

      # Build local operators
      lop[e] = locoperator(SBPp, Nr[e], Ns[e], xt, yt, LFToB = FToB[EToF[:, e]])
    end
    #}}}

    # Assemble the global volume operators
    (M, T, D, vstarts, FToλstarts) =
    LocalGlobalOperators(lop, Nr, Ns, FToB, FToE, FToLF, EToO, EToS,
                         (x) -> cholesky(Symmetric(x)))
    locfactors = M.F

    # Get a unique array indexes for the face to jumps map
    FToδstarts = bcstarts(FToB, FToE, FToLF, (7:1000), Nr, Ns)

    # Compute the number of volume, trace (λ), and jump (δ) points
    VNp = vstarts[nelems+1]-1
    λNp = FToλstarts[nfaces+1]-1
    δNp = FToδstarts[nfaces+1]-1
    # @show (VNp, λNp, δNp)

    # Build the (sparse) λ matrix using the schur complement and factor
    B = assembleλmatrix(FToλstarts, vstarts, EToF, FToB, locfactors, D, T)
    BF = cholesky(Symmetric(B))

    (bλ, λ) = (zeros(λNp), zeros(λNp))
    (Δ, u, g) = (zeros(VNp), zeros(VNp), zeros(VNp))
    δ = zeros(δNp)
    for f = 1:nfaces
      if FToB[f] ∈ (TOP_MAIN, BOTTOM_MAIN, BRANCH)
        (e1, e2) = FToE[:, f]
        (lf1, lf2) = FToLF[:, f]
        (~, ~, L, (x, y), ~, ~, ~, ~, ~, ~, ~) = lop[e1]
        xf = L[lf1] * x
        yf = L[lf1] * y
        @views δ[FToδstarts[f]:(FToδstarts[f+1]-1)] = vex(xf, yf, e2) - vex(xf, yf, e1)
      end
    end

    bc_Dirichlet = (lf, x, y, e, δ) -> vex(x,y,e)
    bc_Neumann   = (lf, x, y, nx, ny, e, δ) -> (nx .* vex_x(x, y, e)
                                                + ny .* vex_y(x, y, e))
    in_jump      = (lf, x, y, e, δ) -> begin
      f = EToF[lf, e]
      if EToS[lf, e] == 1
        if EToO[lf, e]
          return -δ[FToδstarts[f]:(FToδstarts[f+1]-1)]
        else
          return -δ[(FToδstarts[f+1]-1):-1:FToδstarts[f]]
        end
      else
        if EToO[lf, e]
          return  δ[FToδstarts[f]:(FToδstarts[f+1]-1)]
        else
          return  δ[(FToδstarts[f+1]-1):-1:FToδstarts[f]]
        end
      end
    end

    for e = 1:nelems
      locbcarray!((@view g[vstarts[e]:vstarts[e+1]-1]), lop[e], FToB[EToF[:,e]],
                  bc_Dirichlet, bc_Neumann, in_jump, (e, δ))
      (~, ~, ~, (x, y), JH, ~, ~, ~, ~, ~, ~, ~) = lop[e]
      @views g[vstarts[e]:vstarts[e+1]-1] -= JH * fex(x, y, e)

    end

    LocalToGLobalRHS!(bλ, g, u, locfactors, T, vstarts)
    λ[:] = BF \ bλ

    u[:] = T' * λ
    u[:] .= g .+ u

    p1 = plot()
    p2 = plot()
    p3 = plot()
    for e = 1:nelems
      F = locfactors[e]
      (~, ~, ~, (x, y), JH, ~, ~, ~, ~, ~, ~, ~) = lop[e]

      @views u[vstarts[e]:(vstarts[e+1]-1)] = F \ u[vstarts[e]:(vstarts[e+1]-1)]
      #=
      ldiv!((@view u[vstarts[e]:(vstarts[e+1]-1)]), F,
            (@view u[vstarts[e]:(vstarts[e+1]-1)]))
      =#

      @views Δ[vstarts[e]:(vstarts[e+1]-1)] = u[vstarts[e]:(vstarts[e+1]-1)] - vex(x, y, e)
      ϵ[lvl] += Δ[vstarts[e]:(vstarts[e+1]-1)]' * JH * Δ[vstarts[e]:(vstarts[e+1]-1)]
      @plotting begin
        plot!(p1, reshape(x, Nr[e]+1, Ns[e]+1),
              reshape(y, Nr[e]+1, Ns[e]+1),
              reshape(u[vstarts[e]:(vstarts[e+1]-1)], Nr[e]+1, Ns[e]+1))
        plot!(p2, reshape(x, Nr[e]+1, Ns[e]+1),
              reshape(y, Nr[e]+1, Ns[e]+1),
              reshape(vex(x, y, e), Nr[e]+1, Ns[e]+1))
        plot!(p3, reshape(x, Nr[e]+1, Ns[e]+1),
              reshape(y, Nr[e]+1, Ns[e]+1),
              reshape(Δ[vstarts[e]:(vstarts[e+1]-1)], Nr[e]+1, Ns[e]+1))
      end
    end
    ϵ[lvl] = sqrt(ϵ[lvl])
    @show (lvl, ϵ[lvl])
    @plotting begin
      plot!(p1, aspect_ratio = 1, camera = (15, 45), legend=:none)
      plot!(p2, aspect_ratio = 1, camera = (15, 45), legend=:none)
      plot!(p3, aspect_ratio = 1, camera = (15, 45), legend=:none)
      display(plot(p1, p2, p3, layout = (1,3)))
    end
  end
  @show ((log.(ϵ[1:end-1]) - log.(ϵ[2:end])) / log(2))

  nothing
end
