do_plotting = true
if !isdefined(:mtime_global_curved)
  mtime_global_curved = 0
end
if mtime_global_curved < mtime("global_curved.jl")
  println("including global_curved")
  include("global_curved.jl")
  mtime_global_curved = mtime("global_curved.jl")
end
if VERSION <= v"0.6.999999"
  ldiv! = A_ldiv_B!
  cholesky = cholfact
  # using CholmodSolve2
end
using DifferentialEquations

let

  # (verts, EToV, EToF, FToB) = read_inp_2d("meshes/flower_v2.inp")

  #=
  verts = ((-15,   0), (0,   0), (15,   0),
           (-15, -15), (0, -15), (15, -15))
  EToV = ((4, 5, 1, 2), (5, 6, 2, 3))
  EToF = ((1, 2, 3, 4), (2, 5, 6, 7))

  EToV = ((5, 6, 2, 3), (4, 5, 1, 2))
  EToF = ((2, 5, 6, 7), (1, 2, 3, 4))

  EToV = ((6, 3, 5, 2), (2, 1, 5, 4))
  EToF = ((6, 7, 5, 2), (2, 1, 4, 3))

  FToB = (BC_DIRICHLET,
          BC_JUMP_INTERFACE,
          BC_NEUMANN,
          BC_NEUMANN,
          BC_DIRICHLET,
          BC_NEUMANN,
          BC_NEUMANN)
  =#

  verts = ((-15,   0),           (0,   0), (10,  0), (15,   0),
                       (-5, -5), (0,  -5), (10, -5),
           (-15, -15), (-5,-15), (0, -15),           (15, -15))
  EToV = (( 8,  9,  1,  5),
          ( 9, 10,  5,  6),
          (10, 11,  6,  7),
          ( 7, 11,  3,  4),
          ( 6,  7,  2,  3),
          ( 5,  6,  1,  2))
  EToF = ((1, 2, 3, 4),
          (2, 5, 6, 7),
          (5, 8, 9, 10),
          (12, 13, 8, 11),
          (14, 12, 10, 15),
          (4, 14, 7, 16))
  FToB = fill(BC_LOCKED_INTERFACE, (16,))
  for f ∈ (1, 13)
    FToB[f] = BC_DIRICHLET
  end
  for f ∈ (3, 6, 9, 16, 15, 11)
    FToB[f] = BC_NEUMANN
  end
  for f ∈ (5, 14)
    FToB[f] = BC_JUMP_INTERFACE
  end

  #=
  verts = ((-15,   0), (0,   0), (15,   0),
           (-15, -5), (0, -5), (15, -5),
           (-15, -15), (0, -15), (15, -15))
  EToV = ((4, 5, 1, 2),
          (5, 6, 2, 3),
          (7, 8, 4, 5),
          (8, 9, 5, 6))
  EToF = ((1, 2, 3, 4),
          (2, 5, 6, 7),
          (8, 9, 10, 3),
          (9, 11, 12, 6))


  FToB = (BC_DIRICHLET,
          BC_JUMP_INTERFACE,
          BC_LOCKED_INTERFACE,
          BC_NEUMANN,
          BC_DIRICHLET,
          BC_LOCKED_INTERFACE,
          BC_NEUMANN,
          BC_DIRICHLET,
          BC_JUMP_INTERFACE,
          BC_NEUMANN,
          BC_DIRICHLET,
          BC_NEUMANN)
  =#


  verts = flatten_tuples(verts)
  EToV  = flatten_tuples(EToV)
  EToF  = flatten_tuples(EToF)

  # number of elements and faces
  (nelems, nfaces) = (size(EToV, 2), size(FToB, 1))
  println("(nelems, nfaces) = ", (nelems, nfaces))
  #=
  @plotting (p1, p2, p3) = (plot(), plot(), plot())
  @plotting let
    # Do some plotting
    scatter!(p1, verts[1,:], verts[2,:], marker=1, legend=:none)
    display(plot!(p1, aspect_ratio = 1))
    for e = 1:nelems
      V = EToV[:, e]
      plot!(p1, verts[1, V], verts[2, V], linewidth=3)
    end
    display(plot!(p1, aspect_ratio = 1))
  end
  =#

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
    # display(plot!(p1, aspect_ratio = 1))
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

  #{{{ Compute the funtions
  Lx = 15
  hz = 20
  ulinear(x,y) = (hz / Lx) * (x)
  bc_Dirichlet = (lf, x, y, e, δ, t) -> ulinear(x,y)
  bc_Neumann   = (lf, x, y, nx, ny, e, δ, t) -> zeros(size(x))
  in_jump      = (lf, x, y, e, δ, t) -> begin
    f = EToF[lf, e]
    if EToS[lf, e] == 1
      return -δ[FToδstarts[f]:(FToδstarts[f+1]-1)]
    else
      return δ[FToδstarts[f]:(FToδstarts[f+1]-1)]
    end
  end
  #}}}

  #{{{ lithostatic normal stress
  ρ = 3000
  σn = 50*ones(δNp)
  #=
  gravity = 10
  for f = 1:nfaces
    if FToB[f] == BC_JUMP_INTERFACE
      e  = FToE[ 1, f]
      lf = FToLF[1, f]
      (~, ~, L, (~, y), ~, ~, ~, ~, ~, ~, ~) = lop[e]
      δrng = FToδstarts[f]:(FToδstarts[f+1]-1)
      σn[δrng] = -1e-3 * ρ * gravity * (L[lf] * y)
    end
  end
  =#
  μ = 0.6
  shear_modulus = 30
  cs = √(shear_modulus/(ρ * 1e-3))
  η = shear_modulus / (2 * cs)
  #}}}

  # array of jumps
  odefun(V, δ, p, t) = begin
    for e = 1:nelems
      locbcarray!((@view g[vstarts[e]:vstarts[e+1]-1]), lop[e], FToB[EToF[:,e]],
                  bc_Dirichlet, bc_Neumann, in_jump, (e, δ, t))
    end
    LocalToGLobalRHS!(bλ, g, u, locfactors, T, vstarts)
    λ[:] = BF \ bλ

    u[:] = T' * λ
    u[:] .= g .+ u
    for e = 1:nelems
      F = locfactors[e]
      @views u[vstarts[e]:(vstarts[e+1]-1)] = F \ u[vstarts[e]:(vstarts[e+1]-1)]
    end

    # Compute the shear-traction and update velocity
    # mxTz1 = typemin(Float64)
    # mnTz1 = typemax(Float64)
    # mxλ = typemin(Float64)
    # mnλ = typemax(Float64)
    for f = 1:nfaces
      if FToB[f] == BC_JUMP_INTERFACE
        (e1, e2) = FToE[:, f]
        (lf1, lf2) = FToLF[:, f]
        # (~, ~, ~, ~, ~, ~, nx, ~, ~, ~, ~) = lop[e1]
        δrng = FToδstarts[f]:(FToδstarts[f+1]-1)

        Tz1 = shear_modulus * computeTz(f, λ, FToλstarts, u, vstarts, lop, FToE,
                                        FToLF; δ=δ[δrng])
        #=
        Tz2 = computeTz(f, λ, FToλstarts, u, vstarts, lop, FToE, FToLF, EToO;
                        δ=δ[δrng])
        println(Tz1 ≈ -Tz2)
        =#

        # TODO: figure out why this is minus and not plus (error in sign of
        # δ somwhere?)
        V[δrng] = -sign.(Tz1) .* max.(0, (abs.(Tz1) - μ * σn[δrng]) / η)

        # mxTz1 = max(mxTz1, maximum((Tz1)))
        # mnTz1 = min(mnTz1, minimum((Tz1)))
        # λrng = FToλstarts[f]:(FToλstarts[f+1]-1)
        # mxλ = max(mxλ, maximum(abs.(λ[λrng])))
        # mnλ = min(mnλ, minimum(abs.(λ[λrng])))
      end
    end
    # println((t, extrema(V), (mnTz1,mxTz1), (mnλ, mxλ), extrema(abs.(δ))))
    println((t, extrema(V), extrema(abs.(δ))))
    # V[:].=0
    V
  end
  δ = zeros(δNp)

  #=
  for f = 1:nfaces
    if FToB[f] == BC_JUMP_INTERFACE
      (e1, e2) = FToE[:, f]
      (lf1, lf2) = FToLF[:, f]
      (~, ~, ~, ~, ~, ~, nx, ~, ~, ~, ~) = lop[e1]
      δrng = FToδstarts[f]:(FToδstarts[f+1]-1)
      δ[δrng] = -sign.(nx[lf1]) * 90
    end
  end
  =#

  #=
  sol = solve(SteadyStateProblem(odefun, δ))
  δ[:] = sol.u[:]
  =#

  tspan = (0.0, 20.0)
  prob = ODEProblem(odefun, δ, tspan)
  sol = solve(prob,Tsit5(), reltol=1e-8,abstol=1e-8)
  δ = sol.u[end]

  V = zeros(δNp)
  odefun(V, δ, (), tspan[end])
  println("extrema(V) = ", extrema(V))
  println("extrema(δ) = ", extrema(δ))

  @plotting let
    #=
    mx = 0.0
    mn = 0.0
    for e = 1:nelems
      (x, y) = lop[e][4]
      Δu = u[vstarts[e]:(vstarts[e+1]-1)] - ulinear(x,y)
      mn = min(mn, minimum(Δu))
      mx = max(mx, maximum(Δu))
    end
    =#
    mx = maximum(abs.(δ))
    clims = (-mx, mx)
    p2 = plot()
    for e = 1:nelems
      (x, y) = lop[e][4]
      Δu = u[vstarts[e]:(vstarts[e+1]-1)] - ulinear(x,y)
      # Δu = u[vstarts[e]:(vstarts[e+1]-1)]
      plot!(p2, reshape(x, Nr[e]+1, Ns[e]+1),
            reshape(y, Nr[e]+1, Ns[e]+1),
            reshape(Δu, Nr[e]+1, Ns[e]+1),
            st = :surface, c = :balance, clims = clims)
    end
    plot!(p2, aspect_ratio = 1, camera = (15, 26))
    # plot!(p2, aspect_ratio = 1)
    # plot!(p2, aspect_ratio = 1, camera = (45, 45))

    display(plot(p1, p2, layout = (2,1), size = (1400, 800)))
  end
  nothing
end
