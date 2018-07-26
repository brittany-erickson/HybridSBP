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

  (verts, EToV, EToF, FToB) = read_inp_2d("meshes/flower_v2.inp")

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

  #=
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
  =#

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

  if typeof(verts) <: Tuple
    verts = flatten_tuples(verts)
    EToV  = flatten_tuples(EToV)
    EToF  = flatten_tuples(EToF)
  end

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
  hz = 15.1
  ulinear(x, y, t) = (x/Lx) * 1e-9 * t
  bc_Dirichlet = (lf, x, y, e, δ, t) -> (x/Lx) * 1e-9 * t
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

  #{{{ parameters
  ρ = 2.670
  cs = 3.464
  μshear = cs^2 * ρ
  σn = 50 * ones(δNp)
  η = μshear / (2 * cs)
  RSa = 0.025
  RSb = 0.015
  RSDc = 0.008
  RSf0 = 0.6
  RSV0 = 1e-6
  RSVinit = 1e-9
  τz0 = σn .* RSa .* asinh.(RSVinit ./ (2 .* RSV0) .*
                         exp.((RSf0 .+ RSb .* log.(RSV0 ./ RSVinit)) ./ RSa)) .+
            η .* RSVinit

  θ = (RSDc ./ RSV0) .* exp.((RSa ./ RSb) .* log.((2 .* RSV0 ./ RSVinit) .*
      sinh.((τz0 .- η .* RSVinit) ./ (RSa .* σn))) .- RSf0 ./ RSb)
  ψ0 = RSf0 .+ RSb .* log.(RSV0 .* θ ./ RSDc)
  μ = 0.6
  #}}}

  # array of jumps
  odefun(dψV, ψδ, p, t) = begin
    ψ  = @view ψδ[        (1:δNp) ]
    δ  = @view ψδ[ δNp .+ (1:δNp) ]
    dψ = @view dψV[       (1:δNp) ]
    V  = @view dψV[δNp .+ (1:δNp) ]
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
    for f = 1:nfaces
      if FToB[f] == BC_JUMP_INTERFACE
        (e1, e2) = FToE[:, f]
        (lf1, lf2) = FToLF[:, f]
        δrng = FToδstarts[f]:(FToδstarts[f+1]-1)

        Tz = μshear * computeTz(f, λ, FToλstarts, u, vstarts, lop, FToE,
                                        FToLF; δ=δ[δrng])
        for n = 1:length(δrng)
          δn = δrng[n]
          VR = abs(Tz[n] / η)
          VL = -VR
          (Vnew, ~, iter) = newtbndv((V) -> rateandstate(V, ψ[δn], σn[δn],
                                                          Tz[n], η, RSa,
                                                          RSV0),
                                     VL, VR, 0.0)
          if isnan(Vnew) || iter < 0
            println((VL, VR, V[δn], Tz[n], η, RSa, RSV0))
            error()
          end
          if abs(Vnew) > 100
            println(">>>>>>>>>>>", (Vnew, (VL, VR, V[δn], Tz[n], η, RSa, RSV0)))
          end
          V[δn] = Vnew

          dψ(δn) = (RSb * RSV0 / RSDc) * (exp((RSf0-ψ[δn]) / RSb) - abs(V[δn])
                                          / RSV0)
        end

        # TODO: figure out why this is minus and not plus (error in sign of
        # δ somwhere?)
        # V[δrng] = -sign.(Tz) .* max.(0, (abs.(Tz) - μ * σn[δrng]) / η)
      end
    end
    println((t, extrema(V), extrema(abs.(δ))))
    V
  end
  ψδ = zeros(2 * δNp)
  ψδ[1:δNp] .= ψ0

  tspan = (8.8e9, 2.0e10)
  prob = ODEProblem(odefun, ψδ, tspan)
  sol = solve(prob,Tsit5(); dtmax=1e5,dt = 1e-8, atol = 1e-10, rtol = 1e-10, save_everystep=false)
  ψδ = sol.u[end]
  ψδ[δNp .+ (1:δNp)] = ψδ[δNp .+ (1:δNp)]
  δ = @view ψδ[δNp .+ (1:δNp)]

  dψV = zeros(2 * δNp)
  odefun(dψV, ψδ, (), tspan[end])
  V = @view dψV[δNp .+ (1:δNp)]
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
      Δu = u[vstarts[e]:(vstarts[e+1]-1)] - ulinear(x,y,tspan[end])
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
