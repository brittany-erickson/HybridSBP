include("global_curved.jl")

let
  SBPp   = 6 # SBP interior order

  (verts, EToV, EToF, FToB) = read_inp_2d("meshes/flower_v2.inp")

  #=
  Lx = 80
  verts = ((-Lx,   0), (0,   0), (Lx,   0), # 1 2 3
           (-Lx, -40), (0, -40), (Lx, -40), # 4 5 6
           (-Lx, -80), (0, -80), (Lx, -80)) # 7 8 9
  EToV = ((4, 5, 1, 2),
          (5, 6, 2, 3),
          (7, 8, 4, 5),
          (8, 9, 5, 6))
  EToF = ((1,  2,  3, 4),
          (2,  5,  6, 7),
          (8,  9, 10, 3),
          (9, 11, 12, 6))
  FToB = fill(BC_LOCKED_INTERFACE, (12,))
  for f ∈ (1, 5, 8, 11, 4, 7, 10, 12)
    FToB[f] = BC_DIRICHLET
  end
  if typeof(verts) <: Tuple
    verts = flatten_tuples(verts)
    EToV  = flatten_tuples(EToV)
    EToF  = flatten_tuples(EToF)
  end
  =#

  N1 = N0 = 16
  Lx = maximum(verts[1,:])
  Ly = maximum(abs.(verts[2,:]))
  @show (Lx, Ly)

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
  (kx, ky) = (4*π / Lx, 4*π / Ly)
  vex   = (x,y,e) ->       cos.(kx * x) .* cosh.(ky * y)
  vex_x = (x,y,e) -> -kx * sin.(kx * x) .* cosh.(ky * y)
  vex_y = (x,y,e) ->  ky * cos.(kx * x) .* sinh.(ky * y)

  ϵ = zeros(4)
  for lvl = 1:length(ϵ)
    # Set up the local grid dimensions
    Nr = EToN0[1, :] * (2^(lvl-1))
    Ns = EToN0[2, :] * (2^(lvl-1))

    #
    # Build the local volume operators
    #

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

    #
    # Assemble the global volume operators
    #
    (M, FbarT, D, vstarts, FToλstarts) =
    LocalGlobalOperators(lop, Nr, Ns, FToB, FToE, FToLF, EToO, EToS,
                         (x) -> cholesky(Symmetric(x)))
    locfactors = M.F

    # Get a unique array indexes for the face to jumps map
    FToδstarts = bcstarts(FToB, FToE, FToLF, BC_JUMP_INTERFACE, Nr, Ns)

    # Compute the number of volume, trace (λ), and jump (δ) points
    VNp = vstarts[nelems+1]-1
    λNp = FToλstarts[nfaces+1]-1
    δNp = FToδstarts[nfaces+1]-1
    # @show (VNp, λNp, δNp)

    # Build the (sparse) λ matrix using the schur complement and factor
    B = assembleλmatrix(FToλstarts, vstarts, EToF, FToB, locfactors, D, FbarT)
    BF = cholesky(Symmetric(B))

    (bλ, λ, gδ) = (zeros(λNp), zeros(λNp), zeros(λNp))
    (Δ, u, g) = (zeros(VNp), zeros(VNp), zeros(VNp))
    δ = zeros(δNp)

    bc_Dirichlet = (lf, x, y, e, δ) -> vex(x, y, e)
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
      gδe = ntuple(4) do lf
        f = EToF[lf, e]
        if EToO[lf, e]
          return @view gδ[FToλstarts[f]:(FToλstarts[f+1]-1)]
        else
          return  @view gδ[(FToλstarts[f+1]-1):-1:FToλstarts[f]]
        end
      end
      locbcarray!((@view g[vstarts[e]:vstarts[e+1]-1]), gδe, lop[e],
                  FToB[EToF[:,e]], bc_Dirichlet, bc_Neumann, in_jump, (e, δ))
    end

    LocalToGLobalRHS!(bλ, g, gδ,  u, locfactors, FbarT, vstarts)
    #TODO: NEED TO fix for discontinuous τ
    λ[:] = BF \ bλ

    u[:] = -FbarT' * λ
    u[:] .= g .+ u

    for e = 1:nelems
      F = locfactors[e]
      (x, y) = lop[e].coord
      JH = lop[e].JH

      @views u[vstarts[e]:(vstarts[e+1]-1)] = F \ u[vstarts[e]:(vstarts[e+1]-1)]
      #=
      ldiv!((@view u[vstarts[e]:(vstarts[e+1]-1)]), F,
            (@view u[vstarts[e]:(vstarts[e+1]-1)]))
      =#

      @views Δ[vstarts[e]:(vstarts[e+1]-1)] = u[vstarts[e]:(vstarts[e+1]-1)] - vex(x, y, e)
      ϵ[lvl] += Δ[vstarts[e]:(vstarts[e+1]-1)]' * JH * Δ[vstarts[e]:(vstarts[e+1]-1)]
    end
    ϵ[lvl] = sqrt(ϵ[lvl])
    @show (lvl, ϵ[lvl])
  end
  println((log.(ϵ[1:end-1]) - log.(ϵ[2:end])) / log(2))

  nothing
end
