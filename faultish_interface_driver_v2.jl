do_plotting = true
include("global_curved.jl")
if VERSION <= v"0.6.999999"
  ldiv! = A_ldiv_B!
  cholesky = cholfact
  # using CholmodSolve2
end

let

  # This is just needed because functions below expect arrays
  (verts, EToV, EToF, FToB) = read_inp_2d("meshes/BP1_V1.inp")

  # number of elements and faces
  (nelems, nfaces) = (size(EToV, 2), size(FToB, 1))

  # Some sanity checks
  @assert typeof(EToV) == Array{Int, 2} && size(EToV) == (4, nelems)
  @assert typeof(EToF) == Array{Int, 2} && size(EToF) == (4, nelems)
  @assert maximum(maximum(EToF)) == nfaces

  # Determine secondary arrays
  # FToE : Unique Global Face to Element Number
  # FToLF: Unique Global Face to Element local face number
  # EToO : Element to Unique Global Faces Orientation
  # EToS : Element to Unique Global Face Side
  (FToE, FToLF, EToO, EToS) = connectivityarrays(EToV, EToF)

  #{{{ Plot the connectivity using vertices corners
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
      elseif FToB[f] == BC_JUMP_INTERFACE
        plot!(p1, verts[1, V], verts[2, V], color=:green, linewidth=3)
      else
        error("invalid bc")
      end
    end
    plot!(p1, aspect_ratio = 1)
  end
  #}}}

  p   = 4 # SBP interior order
  lvl = 1 # Refinement

  # Set up the local grid dimensions
  EToN0 = fill(13, (2, nelems))
  Nr = EToN0[1, :] * (2^(lvl-1))
  Ns = EToN0[2, :] * (2^(lvl-1))

  #{{{ Build the local volume operators
  # Dictionary to store the operators
  OPTYPE = typeof(locoperator(2, 8, 8, (r,s)->r, (r,s)->s))
  lop = Dict{Int64, OPTYPE}()
  for e = 1:nelems
    # println((e, nelems))
    # Get the element corners
    (x1, x2, x3, x4) = verts[1, EToV[:, e]]
    (y1, y2, y3, y4) = verts[2, EToV[:, e]]

    xt = (r,s)->transfinite_blend(x1, x2, x3, x4, r, s)
    yt = (r,s)->transfinite_blend(y1, y2, y3, y4, r, s)

    # Build local operators
    lop[e] = locoperator(p, Nr[e], Ns[e], xt, yt, LFToB = FToB[EToF[:, e]])
  end
  #}}}

  # Assemble the global volume operators
  (M, T, D, vstarts, FToλstarts) =
  LocalGlobalOperators(lop, Nr, Ns, FToB, FToE, FToLF, EToO, EToS,
                       (x) -> cholesky(Symmetric(x)))
  locfactors = M.F

  # Get a unique array indexes for the face to jumps map
  FToδstarts = bcstarts(FToB, FToE, FToLF, BC_JUMP_INTERFACE, Nr, Ns)

  # Compute the number of volume, trace (λ), and jump (δ) points
  VNp = vstarts[nelems+1]-1
  λNp = FToλstarts[nfaces+1]-1
  δNp = FToδstarts[nfaces+1]-1
  println("(VNp, λNp, δNp) = ", (VNp, λNp, δNp))

  # Build the (sparse) λ matrix using the schur complement and factor
  B = assembleλmatrix(FToλstarts, vstarts, EToF, FToB, locfactors, D, T)
  BF = cholesky(Symmetric(B))

  # Set up some needed arrays
  (bλ, λ, u, g) = (zeros(λNp), zeros(λNp), zeros(VNp), zeros(VNp))

  # array of jumps
  δ = zeros(δNp)

  #{{{ Compute the funtions
  bc_Dirichlet = (lf, x, y, e, t) -> zeros(size(x))
  bc_Neumann   = (lf, x, y, nx, ny, e, t) -> zeros(size(x))
  in_jump      = (lf, x, y, e, t) -> begin
    f = EToF[lf, e]
    if EToS[lf, e] == 1
      return δ[FToδstarts[f]:(FToδstarts[f+1]-1)]
    else
      return -δ[FToδstarts[f]:(FToδstarts[f+1]-1)]
    end
  end
  #}}}

  λ[:] .= 0
  for t = linspace(0,1,10)
  # for t = 0
    δ[:] .= t
    for e = 1:nelems
      locbcarray!((@view g[vstarts[e]:vstarts[e+1]-1]), lop[e], FToB[EToF[:,e]],
                  bc_Dirichlet, bc_Neumann, in_jump, (e,t))
    end
    LocalToGLobalRHS!(bλ, g, u, locfactors, T, vstarts)
    λ[:] = BF \ bλ

    u[:] = T' * λ
    u[:] .= g .+ u
    for e = 1:nelems
      F = locfactors[e]
      @views u[vstarts[e]:(vstarts[e+1]-1)] = F \ u[vstarts[e]:(vstarts[e+1]-1)]
    end

    @plotting let
      clims = (minimum(u), maximum(u))
      p2 = plot()
      for e = 1:nelems
        (x, y) = lop[e][4]
        up = u[vstarts[e]:(vstarts[e+1]-1)]
        plot!(p2, reshape(x, Nr[e]+1, Ns[e]+1),
              reshape(y, Nr[e]+1, Ns[e]+1),
              reshape(up, Nr[e]+1, Ns[e]+1),
              st = :surface, c = :balance, clims = clims)
      end
      # plot!(p2, aspect_ratio = 1, camera = (0, 90))
      # plot!(p2, aspect_ratio = 1)
      plot!(p2, aspect_ratio = 1, camera = (45, 45))

      display(plot(p1, p2, layout = (2,1), size = (2500, 1500)))
      println(t)
      sleep(1)
    end
  end
  nothing
end
