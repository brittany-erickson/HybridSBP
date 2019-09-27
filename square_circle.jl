do_plotting = false
include("global_curved.jl")

let
  RS_FAULT = 7
  VP_FAULT = 8

  SBPp   = 6 # SBP interior order

  # Side set to BC maps
  bc_map = [BC_DIRICHLET, BC_DIRICHLET, BC_NEUMANN, BC_NEUMANN,
            BC_JUMP_INTERFACE]
  (verts, EToV, EToF, FToB, EToBlock) = read_inp_2d("meshes/square_circle.inp";
                                                    bc_map = bc_map)

  for v in 1:size(verts, 2)
    x, y = verts[1, v], verts[2, v]
    if abs(hypot(x,y) - 1) < 1e-5
      Q = atan(y, x)
      verts[1, v], verts[2, v] = cos(Q), sin(Q)
    end

  end
  if typeof(verts) <: Tuple
    verts = flatten_tuples(verts)
    EToV  = flatten_tuples(EToV)
    EToF  = flatten_tuples(EToF)
  end
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
  vex(x,y,e) = begin
    if EToBlock[e] == 1
      return cos.(kx * x) .* cosh.(ky * y)
    elseif EToBlock[e] == 2
      return 10 .+ cos.(kx * x) .* cosh.(ky * y)
    else
      error("invalid block")
    end
  end
  vex_x(x,y,e) = begin
    if EToBlock[e] == 1
      return -kx * sin.(kx * x) .* cosh.(ky * y)
    elseif EToBlock[e] == 2
      return -kx * sin.(kx * x) .* cosh.(ky * y)
    else
      error("invalid block")
    end
  end
  vex_y(x,y,e) = begin
    if EToBlock[e] == 1
      return ky * cos.(kx * x) .* sinh.(ky * y)
    elseif EToBlock[e] == 2
      return ky * cos.(kx * x) .* sinh.(ky * y)
    else
      error("invalid block")
    end
  end

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

      ex = [(α) -> x1 * (1 .- α) / 2 + x3 * (1 .+ α) / 2,
            (α) -> x2 * (1 .- α) / 2 + x4 * (1 .+ α) / 2,
            (α) -> x1 * (1 .- α) / 2 + x2 * (1 .+ α) / 2,
            (α) -> x3 * (1 .- α) / 2 + x4 * (1 .+ α) / 2]
      ey = [(α) -> y1 * (1 .- α) / 2 + y3 * (1 .+ α) / 2,
            (α) -> y2 * (1 .- α) / 2 + y4 * (1 .+ α) / 2,
            (α) -> y1 * (1 .- α) / 2 + y2 * (1 .+ α) / 2,
            (α) -> y3 * (1 .- α) / 2 + y4 * (1 .+ α) / 2]
      if FToB[EToF[1, e]] == BC_JUMP_INTERFACE
          error("fix me")
      end
      if FToB[EToF[2, e]] == BC_JUMP_INTERFACE
          error("fix me")
      end
      if FToB[EToF[3, e]] == BC_JUMP_INTERFACE
        Q1 = atan(y1, x1)
        Q2 = atan(y2, x2)
        ex[3] = (α) -> cos.(Q1 * (1 .- α) / 2 + Q2 * (1 .+ α) / 2)
        ey[3] = (α) -> sin.(Q1 * (1 .- α) / 2 + Q2 * (1 .+ α) / 2)
        if !(-π/2 < Q1 - Q2 < π/2)
          Q2 -= sign(Q2) * 2 * π
        end
      end
      if FToB[EToF[4, e]] == BC_JUMP_INTERFACE
        Q3 = atan(y3, x3)
        Q4 = atan(y4, x4)
        ex[4] = (α) -> cos.(Q3 * (1 .- α) / 2 + Q4 * (1 .+ α) / 2)
        ey[4] = (α) -> sin.(Q3 * (1 .- α) / 2 + Q4 * (1 .+ α) / 2)
        if !(-π/2 < Q3 - Q4 < π/2)
          error("fix me")
        end
      end
      xt = (r,s)->transfinite_blend(x1, x2, x3, x4, r, s;
                                    e1=ex[1], e2=ex[2], e3=ex[3], e4=ex[4])
      yt = (r,s)->transfinite_blend(y1, y2, y3, y4, r, s;
                                    e1=ey[1], e2=ey[2], e3=ey[3], e4=ey[4])

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
    FToδstarts = bcstarts(FToB, FToE, FToLF, (RS_FAULT,VP_FAULT), Nr, Ns)

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
    for f = 1:nfaces
      if FToB[f] == BC_JUMP_INTERFACE
        (e1, e2) = FToE[:, f]
        (lf1, lf2) = FToLF[:, f]
        (x, y) = lop[e1].coord
        L = lop[e1].L
        xf = L[lf1] * x
        yf = L[lf1] * y
        @views δ[FToδstarts[f]:(FToδstarts[f+1]-1)] =
            vex(xf, yf, e2) - vex(xf, yf, e1)
      end
    end

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

    lockedblock = Array{Bool, 1}(undef, nelems)
    for e = 1:nelems
      if @views FToB[EToF[:,e]] == [BC_LOCKED_INTERFACE, BC_LOCKED_INTERFACE,
                                    BC_LOCKED_INTERFACE, BC_LOCKED_INTERFACE]
        lockedblock[e] = true
      else
        lockedblock[e] = false
      end
    end
    LocalToGLobalRHS!(bλ, g, gδ,  u, locfactors, FbarT, vstarts, lockedblock)
    #TODO: NEED TO fix for discontinuous τ
    λ[:] = BF \ bλ

    u[:] = -FbarT' * λ
    u[:] .= g .+ u

    @plotting begin
      p1 = plot()
      p2 = plot()
      p3 = plot()
    end
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
      #=
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
      =#
    end
    #=
    @plotting begin
      plot!(p1, aspect_ratio = 1, camera = (15, 26), legend=:none)
      plot!(p2, aspect_ratio = 1, camera = (15, 26), legend=:none)
      plot!(p3, aspect_ratio = 1, camera = (15, 26), legend=:none)
      display(plot(p1, p2, p3, layout = (1,3)))
    end
    =#
    ϵ[lvl] = sqrt(ϵ[lvl])
    @show (lvl, ϵ[lvl])
  end
  println((log.(ϵ[1:end-1]) - log.(ϵ[2:end])) / log(2))

  nothing
end
