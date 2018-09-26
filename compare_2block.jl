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

using Compat.Printf: @sprintf

let
  SBPp    = 2 # SBP interior order
  τscale  = 12

  verts = ((-3, 3), (0, 3), (3, 3),  # 1 2 3
           (-3, 0), (0, 0), (3, 0))  # 4 5 6
  EToV = ((4, 5, 1, 2),
          (5, 6, 2, 3))
  EToF = ((1,  2,  3, 4),
          (2,  5,  6, 7))
  FToB = fill(BC_LOCKED_INTERFACE, (7,))
  EToBlock = (1,2)
  for f ∈ (1, 3, 4, 5, 6, 7)
    FToB[f] = BC_DIRICHLET
  end
  for f ∈ (2,)
    FToB[f] = BC_JUMP_INTERFACE
  end

  EToV = ((5, 2, 4, 1),
          (6, 3, 5, 2))
  EToF = ((1,  2,  3, 4),
          (5,  6,  7, 3))
  FToB = fill(BC_LOCKED_INTERFACE, (7,))
  EToBlock = (1,2)
  for f ∈ (1, 2, 4, 5, 6, 7)
    FToB[f] = BC_DIRICHLET
  end
  for f ∈ (3,)
    FToB[f] = BC_JUMP_INTERFACE
  end
  for f ∈ (1,2,5,6,)
    FToB[f] = BC_NEUMANN
  end

  N0 = 2
  N1 = 2
  lvl = 1 # Refinement


  if typeof(verts) <: Tuple
    verts = flatten_tuples(verts)
    EToV  = flatten_tuples(EToV)
    EToF  = flatten_tuples(EToF)
    EToBlock  = flatten_tuples(EToBlock)
  end

  # number of elements and faces
  (nelems, nfaces) = (size(EToV, 2), size(FToB, 1))
  @show (nelems, nfaces)
  FToB2 = copy(FToB)

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
  @plotting (p1, p2, p3) = (plot(), plot(), plot())
  @plotting let
    # Do some plotting
    scatter!(p1, verts[1,:], verts[2,:], marker=1, legend=:none)
    LFToLV = flatten_tuples(((1,3), (2, 4), (1,2), (3,4)))
    for f = 1:nfaces
      e  = FToE[1,f]
      lf = FToLF[1,f]
      V = EToV[LFToLV[:,lf], e]
      if FToB[f] == BC_DIRICHLET
        plot!(p1, verts[1, V], verts[2, V], color=:red, linewidth=3)
      elseif FToB[f] == BC_NEUMANN
        plot!(p1, verts[1, V], verts[2, V], color=:blue, linewidth=3)
      elseif FToB[f] == BC_LOCKED_INTERFACE
        plot!(p1, verts[1, V], verts[2, V], color=:black, linewidth=1)
      elseif FToB[f] >= BC_JUMP_INTERFACE
        plot!(p1, verts[1, V], verts[2, V], color=:green, linewidth=3)
      else
        error("invalid bc")
      end
    end
    plot!(p1, aspect_ratio = 1)
    display(plot!(p1, aspect_ratio = 1))
  end

  begin
    # Set up the local grid dimensions
    Nr = EToN0[1, :] * (2^(lvl-1))
    Ns = EToN0[2, :] * (2^(lvl-1))

    #{{{ Build the local volume operators
    # Dictionary to store the operators
    OPTYPE = typeof(locoperator(2, 8, 8, (r,s)->r, (r,s)->s; τscale=τscale))
    lop = Dict{Int64, OPTYPE}()
    for e = 1:nelems
      # @show (e, nelems)
      # Get the element corners
      (x1, x2, x3, x4) = verts[1, EToV[:, e]]
      (y1, y2, y3, y4) = verts[2, EToV[:, e]]

      xt = (r,s)->transfinite_blend(x1, x2, x3, x4, r, s)
      yt = (r,s)->transfinite_blend(y1, y2, y3, y4, r, s)

      # Build local operators
      lop[e] = locoperator(SBPp, Nr[e], Ns[e], xt, yt, LFToB = FToB[EToF[:, e]];
                           τscale=τscale, pm = SBPp)
      println("Local Operator for block $(e)")
      display(full(lop[e][1]))
      println()
    end
    #}}}

    # Assemble the global volume operators
    (M, T, D, vstarts, FToλstarts) =
    LocalGlobalOperators(lop, Nr, Ns, FToB, FToE, FToLF, EToO, EToS,
                         (x) -> cholesky(Symmetric(x)))
                         # (x) -> lufact(x))
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
    println("Schur complement matrix")
    display(full(B))
    return
    BF = cholesky(Symmetric(B))

    (bλ, λ) = (zeros(λNp), zeros(λNp))
    (Δ, u, g) = (zeros(VNp), zeros(VNp), zeros(VNp))
    δ = zeros(δNp)
    for f = 1:nfaces
      if FToB[f] ∈ (BC_JUMP_INTERFACE,)
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
    p4 = plot(legend=:none, title="T")
    p5 = plot(legend=:none, title="dT")
    p6 = plot(legend=:none, title="u")
    p7 = plot(legend=:none, title="du")
    lw = 6
    for f = 1:nfaces
      if FToB2[f] ∈ (BC_JUMP_INTERFACE,)
        (e1, e2) = FToE[:, f]
        (lf1, lf2) = FToLF[:, f]
        (~, ~, L1, (x, y), ~, ~, nx, ny, ~, ~, ~) = lop[e1]
        xf = L1[lf1] * x
        yf = L1[lf1] * y
        (~, ~, L2, (x, y), ~, ~, nx, ny, ~, ~, ~) = lop[e2]

        μ = 32
        τex = μ * vex_x(xf,yf,e1) .* nx[lf1] + vex_y(xf,yf,e1) .* ny[lf1]
        (τm, τp) = computeTz4(f, λ, FToλstarts, u, vstarts, lop, FToE, FToLF,
                              EToO)
        τ = μ * (τm .+ τp) / 2
        #=
        δrng = FToδstarts[f]:(FToδstarts[f+1]-1)
        τ = computeTz1(f, λ, FToλstarts, u, vstarts, lop, FToE, FToLF; δ=δ[δrng])
        =#
        # τex = computeTz1(f, λ, FToλstarts, u, vstarts, lop, FToE, FToLF, EToO; δ=δ[δrng])
        # τex = computeTz(f, λ, FToλstarts, u, vstarts, lop, FToE, FToLF; δ=δ[δrng])
        plot!(p4, yf, τex, linewidth = lw)
        plot!(p4, yf, τ, marker=5, linewidth = 0)
        plot!(p5, yf, abs.(τex-τ), linewidth = lw)

        wex = vex(xf,yf,e1) - vex(xf,yf,e2)
        vu1 = @view u[vstarts[e1]:(vstarts[e1+1]-1)]
        vu2 = @view u[vstarts[e2]:(vstarts[e2+1]-1)]
        w = (L1[lf1] * vu1) - (L2[lf2] * vu2)
        plot!(p6, yf, wex, linewidth = lw)
        plot!(p6, yf, w, marker= 5, linewidth = 0)
        plot!(p7, yf, abs.(wex-w), linewidth = lw)
      end
    end
    @plotting begin
      plot!(p1, aspect_ratio = 1, camera = (15, 45), legend=:none, title="u")
      plot!(p2, aspect_ratio = 1, camera = (15, 45), legend=:none, title="uex")
      plot!(p3, aspect_ratio = 1, camera = (15, 45), legend=:none, title="Δu")
      # display(plot(p1, p2, p3, layout = (1,3)))

      display(plot(p4, p6, p5, p7, layout = (2,2), size = (1400, 1400)))
      savefig((@sprintf "compare_p%d_lvl%d_pt%e.pdf" SBPp lvl τscale))
    end
  end
  @show ((log.(ϵ[1:end-1]) - log.(ϵ[2:end])) / log(2))

  nothing
end
