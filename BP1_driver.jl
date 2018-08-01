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
end
using DifferentialEquations
macro expr_println(ex)
  return :(println($(string(ex)), " = ", $(esc(ex))))
end

let
  sim_years = 1348.369
  year_seconds = 31556926.

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
  for f ∈ (1, 5, 8, 11)
    FToB[f] = BC_DIRICHLET
  end
  for f ∈ (4, 7, 10, 12)
    FToB[f] = BC_NEUMANN
  end
  for f ∈ (2, 9)
    FToB[f] = BC_JUMP_INTERFACE
  end

  if typeof(verts) <: Tuple
    verts = flatten_tuples(verts)
    EToV  = flatten_tuples(EToV)
    EToF  = flatten_tuples(EToF)
  end

  # (Nx, Ny1, Ny2) = (90, 40, 41)
  (Nx, Ny1, Ny2) = (13, 30, 20)
  EToN0 = flatten_tuples(((Nx, Ny1), (Nx, Ny1), (Nx, Ny2), (Nx, Ny2)))

  # number of elements and faces
  (nelems, nfaces) = (size(EToV, 2), size(FToB, 1))
  @expr_println (nelems, nfaces)

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
    display(plot!(p1, aspect_ratio = 1))
  end
  #}}}

  p   = 4 # SBP interior order
  lvl = 1 # Refinement

  # Set up the local grid dimensions
  Nr = EToN0[1, :] * (2^(lvl-1))
  Ns = EToN0[2, :] * (2^(lvl-1))

  #{{{ Build the local volume operators
  # Dictionary to store the operators
  OPTYPE = typeof(locoperator(2, 8, 8, (r,s)->r, (r,s)->s))
  lop = Dict{Int64, OPTYPE}()
  for e = 1:nelems
    # @expr_println (e, nelems)
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
  @expr_println (VNp, λNp, δNp)

  # Build the (sparse) λ matrix using the schur complement and factor
  B = assembleλmatrix(FToλstarts, vstarts, EToF, FToB, locfactors, D, T)
  BF = cholesky(Symmetric(B))

  # Set up some needed arrays
  (bλ, λ, u, g) = (zeros(λNp), zeros(λNp), zeros(VNp), zeros(VNp))

  #{{{ Compute the boundary/interface functions
  Vp = 1e-9
  ulinear(x, y, t) = (x/Lx) * (Vp/2) * t
  bc_Dirichlet = (lf, x, y, e, δ, t) -> ulinear(x,y,t)
  bc_Neumann   = (lf, x, y, nx, ny, e, δ, t) -> zeros(size(x))
  in_jump      = (lf, x, y, e, δ, t) -> begin
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
  #}}}

  #{{{ parameters
  ρ = 2.670
  cs = 3.464
  μshear = cs^2 * ρ
  σn = 50 * ones(δNp)
  η = μshear / (2 * cs)
  RSa = 0.010
  RSamax = 0.025
  RSb = 0.015
  RSDc = 0.008
  RSf0 = 0.6
  RSV0 = 1e-6
  RSVinit = 1e-9
  τz0 = σn .* RSamax .* asinh.(RSVinit ./ (2 .* RSV0) .*
                         exp.((RSf0 .+ RSb .* log.(RSV0 ./ RSVinit)) ./ RSamax)) .+
            η .* RSVinit

  θ = (RSDc ./ RSV0) .* exp.((RSa ./ RSb) .* log.((2 .* RSV0 ./ RSVinit) .*
      sinh.((τz0 .- η .* RSVinit) ./ (RSa .* σn))) .- RSf0 ./ RSb)
  ψ0 = RSf0 .+ RSb .* log.(RSV0 .* θ ./ RSDc)

  for f = 1:nfaces
    if FToB[f] == BC_JUMP_INTERFACE
      (e1, e2) = FToE[:, f]
      (lf1, lf2) = FToLF[:, f]
      (~, ~, ~, ~, ~, ~, nx, ~, ~, ~, ~) = lop[e1]
      δrng = FToδstarts[f]:(FToδstarts[f+1]-1)
      for n = 1:length(δrng)
        τz0[δrng[n]] = sign(nx[lf1][n])*abs(τz0[δrng[n]])
      end
    end
  end
  #}}}

  reject_step = [false]
  # array of jumps
  odefun(dψV, ψδ, p, t) = begin
    if reject_step[1]
      return
    end
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
      #=
      ldiv!((@view u[vstarts[e]:(vstarts[e+1]-1)]), F,
            (@view u[vstarts[e]:(vstarts[e+1]-1)]))
            =#
    end

    # Compute the shear-traction and update velocity
    show_val = true
    for f = 1:nfaces
      if FToB[f] == BC_JUMP_INTERFACE
        (e1, e2) = FToE[:, f]
        (lf1, lf2) = FToLF[:, f]
        δrng = FToδstarts[f]:(FToδstarts[f+1]-1)

        #=
        Tz = μshear * computeTz(f, λ, FToλstarts, u, vstarts, lop, FToE,
                                 FToLF; δ=δ[δrng])
        =#
        Tz = μshear * computeTz3(f, λ, FToλstarts, u, vstarts, lop, FToE,
                                 FToLF, EToO)
        for n = 1:length(δrng)
          δn = δrng[n]
          τ = Tz[n] + τz0[δn]
          VR = abs(τ / η)
          VL = -VR
          (Vnew, ~, iter) = newtbndv((V) -> rateandstate(V, ψ[δn], σn[δn],
                                                         τ, η, RSa, RSV0),
                                     VL, VR, 0.0)
          if show_val
            show_val = false
            @expr_println (ψ[δn], σn[δn], τ, η, RSa, RSV0)
          end
          if isnan(Vnew) || iter < 0
            @expr_println (VL, VR, V[δn], Vnew, Tz[n], η, RSa, RSV0)
            Vnew = 1e10
            reject_step[1] = true
            return
            #error()
          end
          #=
          =#
          #=
          if abs(Vnew) > 100
            @expr_println (Vnew, (VL, VR, V[δn], Tz[n], η, RSa, RSV0))
          end
          =#
          V[δn] = Vnew

          dψ[δn] = (RSb * RSV0 / RSDc) * (exp((RSf0-ψ[δn]) / RSb) - abs(V[δn])
                                          / RSV0)
        end
      end
    end
    @expr_println (t/year_seconds, extrema(V), extrema(abs.(δ)))
    V
  end
  ψδ = zeros(2 * δNp)
  ψδ[1:δNp] .= ψ0

  tspan = (0, sim_years * year_seconds)
  prob = ODEProblem(odefun, ψδ, tspan)
  Vmin = 0.0
  δmax = -1
  δcheck(ψδ, p, t) = begin
    δ = @view ψδ[δNp .+ (1:δNp)]
    δnewmax = maximum(abs.(δ))
    if δmax < 0
      δmax = δnewmax
    end
    if reject_step[1]
      reject_step[1] = false
      println(">>>>>>>>>>REJECT (Bool)<<<<<<<<<")
      return true
    end
    if isnan(δnewmax)
      println(">>>>>>>>>>REJECT (NaN)<<<<<<<<<")
      return true
    end
    if δnewmax > 1e-5 && 2*δmax < δnewmax
      println(">>>>>>>>>>REJECT<<<<<<<<<")
      return true
    end
    δmax = δnewmax
    println(">>>>>>>>>>PASS<<<<<<<<<")
    return false
  end
  sol = solve(prob,Tsit5(); isoutofdomain=δcheck, dt = 1e3,
              atol = 1e-10, rtol = 1e-10, save_everystep=false)
  ψδ = sol.u[end]
  ψδ[δNp .+ (1:δNp)] = ψδ[δNp .+ (1:δNp)]
  δ = @view ψδ[δNp .+ (1:δNp)]

  dψV = zeros(2 * δNp)
  odefun(dψV, ψδ, (), tspan[end])
  V = @view dψV[δNp .+ (1:δNp)]
  @expr_println extrema(V)
  @expr_println extrema(δ)

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
