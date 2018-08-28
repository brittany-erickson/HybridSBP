do_plotting = false
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

using DifferentialEquations
using Compat.Printf: @sprintf

let
  RS_FAULT = 7
  VP_FAULT = 8
  sim_years = 3000.
  year_seconds = 31556926.
  SBPp   = 4 # SBP interior order

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
  for f ∈ (2,)
    FToB[f] = RS_FAULT
  end
  for f ∈ (9,)
    FToB[f] = VP_FAULT
  end
  N0 = 1600
  N1 = 800
  lvl = 1 # Refinement
  base_name = "BP1_uniform"

  #=
  (verts, EToV, EToF, FToB) = read_inp_2d("meshes/BP1_V1.inp")
  Lx = maximum(verts[1,:])
  N1 = N0 = 13
  lvl = 3 # Refinement
  base_name = "BP1_mesh_$(lvl)"
  base_name = "BP1_V1_$(lvl)"
  =#

  #=
  (verts, EToV, EToF, FToB) = read_inp_2d("meshes/BP1_V0.inp")
  Lx = maximum(verts[1,:])
  Ly = maximum(abs.(verts[2,:]))
  N1 = N0 = 50
  lvl = 1 # Refinement
  base_name = "BP1_V0_p_$(SBPp)_lvl_$(lvl)"

  #=
  r = verts[1,:]
  s = verts[2,:]
  x = @view verts[1,:]
  y = @view verts[2,:]
  x .= x .+ 3 * sin.(2 * π * s / Ly) .* sin.(2 * π * r / Lx)
  y .= y .+ 3 * sin.(5*π * s / Ly) .* sin.(2 * π * r / Lx)
  base_name = "BP1_V0_skew_p_$(SBPp)_lvl_$(lvl)"
  =#
  =#

  #=
  (verts, EToV, EToF, FToB) = read_inp_2d("meshes/BP1_V2.inp")
  Lx = maximum(verts[1,:])
  N1 = N0 = 50
  lvl = 1 # Refinement
  base_name = "BP1_V2_$(lvl)"
  =#

  if typeof(verts) <: Tuple
    verts = flatten_tuples(verts)
    EToV  = flatten_tuples(EToV)
    EToF  = flatten_tuples(EToF)
  end

  # number of elements and faces
  (nelems, nfaces) = (size(EToV, 2), size(FToB, 1))
  @show (nelems, nfaces)

  EToN0 = zeros(Int64, 2, nelems)
  EToN0[1, :] .= N0
  EToN0[2, :] .= N1

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
      elseif FToB[f] >= BC_JUMP_INTERFACE
        plot!(p1, verts[1, V], verts[2, V], color=:green, linewidth=3)
      else
        error("invalid bc")
      end
    end
    plot!(p1, aspect_ratio = 1)
    display(plot!(p1, aspect_ratio = 1))
  end
  #}}}


  # Set up the local grid dimensions
  Nr = EToN0[1, :] * (2^(lvl-1))
  Ns = EToN0[2, :] * (2^(lvl-1))

  #{{{ Build the local volume operators
  # Dictionary to store the operators
  OPTYPE = typeof(locoperator(2, 8, 8, (r,s)->r, (r,s)->s))
  lop = Dict{Int64, OPTYPE}()
  for e = 1:nelems
    @show (e, nelems)
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
  (vstarts, M, H, X, Y, E) = glovoloperator(lop, Nr, Ns)

  # Build the trace operators
  (FToλstarts, T, D) = gloλoperator(lop, vstarts, FToB, FToE, FToLF, EToO, EToS,
                                 Nr, Ns)
  λNp = FToλstarts[nfaces+1]-1

  # Get a unique array indexes for the face to jumps map
  FToδstarts = bcstarts(FToB, FToE, FToLF, (RS_FAULT,VP_FAULT), Nr, Ns)

  # Compute the number of volume, trace (λ), and jump (δ) points
  VNp = vstarts[nelems+1]-1
  λNp = FToλstarts[nfaces+1]-1
  δNp = FToδstarts[nfaces+1]-1
  @show (VNp, λNp, δNp)

  A = [ M -T'; -T  Diagonal(D) ]
  AF = cholesky(Symmetric(A))

  # Set up some needed arrays
  (uλ, g) = (zeros(VNp + λNp), zeros(VNp+λNp))
  u = @view uλ[1:VNp]
  λ = @view uλ[VNp+(1:λNp)]

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
  RSamin = 0.010
  RSamax = 0.025
  RSb = 0.015
  RSDc = 0.008
  RSf0 = 0.6
  RSV0 = 1e-6
  RSVinit = 1e-9
  RSa = zeros(δNp)
  RSH1 = 15;
  RSH2 = 18;
  fault_y = zeros(δNp)
  for f = 1:nfaces
    if FToB[f] ∈ (RS_FAULT, VP_FAULT)
      (e1, ~) = FToE[:, f]
      (lf1, ~) = FToLF[:, f]
      (~, ~, L, (~, y), ~, ~, ~, ~, ~, ~, ~) = lop[e1]
      fault_y[FToδstarts[f]:(FToδstarts[f+1]-1)] = yf = L[lf1] * y
      δrng = FToδstarts[f]:(FToδstarts[f+1]-1)
      for n = 1:length(δrng)
        RSa[δrng[n]] = RSamin - (RSamin - RSamax) *
          min(1, max(0, (RSH1 + yf[n])/(RSH1 - RSH2)))
      end
    end
  end
  τz0 = σn .* RSamax .* asinh.(RSVinit ./ (2 .* RSV0) .*
                               exp.((RSf0 .+ RSb .* log.(RSV0 ./ RSVinit)) ./
                                    RSamax)) .+ η .* RSVinit

  θ = (RSDc ./ RSV0) .* exp.((RSa ./ RSb) .* log.((2 .* RSV0 ./ RSVinit) .*
      sinh.((τz0 .- η .* RSVinit) ./ (RSa .* σn))) .- RSf0 ./ RSb)
  ψ0 = RSf0 .+ RSb .* log.(RSV0 .* θ ./ RSDc)

  for f = 1:nfaces
    if FToB[f] == RS_FAULT
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

  #{{{ ODE fun
  reject_step = [false]
  # array of jumps
  lockedblock = Array{Bool, 1}(undef, nelems)
  for e = 1:nelems
    if @views FToB[EToF[:,e]] == [BC_LOCKED_INTERFACE, BC_LOCKED_INTERFACE,
                                  BC_LOCKED_INTERFACE, BC_LOCKED_INTERFACE]
      lockedblock[e] = true
    else
      lockedblock[e] = false
    end
  end
  τ = zeros(δNp)
  τp = zeros(δNp)
  τm = zeros(δNp)
  odefun(dψV, ψδ, p, t) = begin
    begin
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
      uλ[:] = AF \ g

      # Compute the shear-traction and update velocity
      show_val = false
      for f = 1:nfaces
        if FToB[f] == RS_FAULT
          (e1, e2) = FToE[:, f]
          (lf1, lf2) = FToLF[:, f]
          δrng = FToδstarts[f]:(FToδstarts[f+1]-1)

          (τm[δrng], τp[δrng]) = computeTz5(f, FToλstarts, u, vstarts, lop,
                                            FToE, FToLF, EToO)
          (τm[δrng], τp[δrng]) = (μshear * τm[δrng], μshear * τp[δrng])
          τ[δrng] .= (τm[δrng] .+ τp[δrng]) / 2

          for n = 1:length(δrng)
            δn = δrng[n]
            τ[δn] .= τ[δn] .+ τz0[δn]
            if isnan(τ[δn])
              reject_step[1] = true
              return
            end
            VR = abs(τ[δn] / η)
            VL = -VR
            (Vnew, ~, iter) = newtbndv((V) -> rateandstate(V, ψ[δn], σn[δn],
                                                           τ[δn], η, RSa[δn],
                                                           RSV0),
                                       VL, VR, 0.0; atolx=1e-9, rtolx=1e-9,
                                       ftol=1e-9)
            if show_val
              show_val = false
              @show (ψ[δn], σn[δn], τ[δn], η, RSa[δn], RSV0)
            end
            if isnan(Vnew) || iter < 0
              # @show (VL, VR, V[δn], Vnew, Tz[n], η, RSa[δn], RSV0)
              println("V reject")
              Vnew = 1e10
              reject_step[1] = true
              return
              #error()
            end
            #=
            =#
            #=
            if abs(Vnew) > 100
              @show (Vnew, (VL, VR, V[δn], Tz[n], η, RSa[δn], RSV0))
            end
            =#
            V[δn] = Vnew

            dψ[δn] = (RSb * RSV0 / RSDc) * (exp((RSf0-ψ[δn]) / RSb) - abs(V[δn])
                                            / RSV0)
            if !isfinite(dψ[δn])
              println("ψ reject")
              dψ[δn] = 0
              reject_step[1] = true
              return
            end

          end
        elseif FToB[f] == VP_FAULT
          (e1, ~) = FToE[:, f]
          (lf1, ~) = FToLF[:, f]
          (~, ~, ~, ~, ~, ~, nx, ~, ~, ~, ~) = lop[e1]
          for δn = FToδstarts[f]:(FToδstarts[f+1]-1)
            V[δn] = sign(nx[lf1][1]) * Vp
          end
        end
      end
      # @show (t/year_seconds, extrema(V), extrema(abs.(δ)), extrema(abs.(ψ)), extrema(dψ))
      V
    end
  end
  ψδ = zeros(2 * δNp)
  ψδ[1:δNp] .= ψ0
  dψV = zeros(2 * δNp)
  #}}}

  #{{{ Setup station data
  stations = [0, -2.5, -5, -7.5, -10, -12.5, -15, -17.5, -20, -25, -30]
  numstations = length(stations)
  station_ind = zeros(Int64, numstations)
  station_face = zeros(Int64, numstations)
  for s = 1:numstations
    n = station_ind[s] = argmin(abs.(fault_y-stations[s]))
    f = station_face[s] = findlast((m) -> m <= n, FToδstarts)
    println((stations[s], fault_y[station_ind[s]], FToδstarts[f] <= n < FToδstarts[f+1]))
  end
  station_t = zeros(Float64, ceil(Int64, sim_years)*10)
  station_V = zeros(Float64, numstations, ceil(Int64, sim_years)*10)
  station_τ = zeros(Float64, numstations, ceil(Int64, sim_years)*10)
  station_θ = zeros(Float64, numstations, ceil(Int64, sim_years)*10)
  station_δ = zeros(Float64, numstations, ceil(Int64, sim_years)*10)
  #}}}

  tspan = (0, sim_years * year_seconds)

  dψV = zeros(2 * δNp)
  # @time odefun(dψV, ψδ, (), 1.0)
  # @time odefun(dψV, ψδ, (), 1.0)
  # @time odefun(dψV, ψδ, (), 1.0)
  # @time odefun(dψV, ψδ, (), 1.0)
  # return
  prob = ODEProblem(odefun, ψδ, tspan)
  Vmin = 0.0
  δmax = -1
  stepcheck(ψδ, p, t) = begin
    if reject_step[1]
      reject_step[1] = false
      println("reject")
      return true
    end
    return false
  end
  tlast = 0
  tdump = 100
  nxt_ind = 1
  fault_order = sortperm(fault_y)
  open("$(base_name)_slip.dat", "w") do f
    write(f, "tsec                   tyear                  ")
    for k = 1:length(fault_order)
      @printf f "%+.16e " fault_y[fault_order[k]]
    end
    write(f, "\n")
  end
  cb = SavingCallback((ψδ, t, i) -> begin
                        Vmax = 0.0
                        if isdefined(i, :fsallast)
                          dψV = i.fsallast
                          V  = @view dψV[δNp .+ (1:δNp) ]
                          Vmax = maximum(abs.(extrema(V)))
                          tnext = tlast + (Vmax > 1e-3 ? 0.1 : year_seconds)
                          if (t >= tnext)
                            ψ  = @view ψδ[        (1:δNp) ]
                            δ  = @view ψδ[ δNp .+ (1:δNp) ]
                            dψ = @view dψV[       (1:δNp) ]
                            tlast = tnext
                            @show (t/year_seconds, Vmax)
                            @show base_name
                            station_t[nxt_ind] = t
                              open("$(base_name)_slip.dat", "a") do f
                                @printf f "%.16e %.16e " t t/year_seconds
                                for k = 1:length(fault_order)
                                  @printf f "%+.16e " δ[fault_order[k]]
                                end
                                write(f, "\n")
                              end
                            for s = 1:numstations
                              n = station_ind[s]
                              station_V[s, nxt_ind] = V[n]
                              station_θ[s, nxt_ind] = RSDc * exp((ψ[n] - RSf0) / RSb) / RSV0
                              station_δ[s, nxt_ind] = δ[n]
                              station_τ[s, nxt_ind] = τ[n] - η * V[n]
                            end
                            @plotting let
                              p1 = plot(fault_y, τ)
                              p2 = plot(fault_y, V)
                              p3 = plot(fault_y, δ)
                              p4 = plot(fault_y, τp - τm)
                              display(plot(p1, p2, p3, p4, layout = (4,1),
                                           size = (1400, 1400),
                                           title = "time = $(t/year_seconds)"))
                              savefig((@sprintf "figs/%s_%08d.png" base_name nxt_ind))
                            end
                            nxt_ind = nxt_ind + 1
                            if nxt_ind > length(station_t)
                              old_len = length(station_t)
                              len = 2 * old_len
                              old_station_V = station_V
                              old_station_τ = station_τ
                              old_station_ψ = station_θ
                              old_station_δ = station_δ
                              station_t = Array{Float64, 1}(undef, len)
                              station_V = Array{Float64, 2}(undef, numstations,  len)
                              station_τ = Array{Float64, 2}(undef, numstations,  len)
                              station_θ = Array{Float64, 2}(undef, numstations,  len)
                              station_δ = Array{Float64, 2}(undef, numstations,  len)
                              station_t[1:old_len] .= old_station_t[1:old_len]
                              station_V[:, 1:old_len] .= old_station_V[:, 1:old_len]
                              station_τ[:, 1:old_len] .= old_station_τ[:, 1:old_len]
                              station_θ[:, 1:old_len] .= old_station_ψ[:, 1:old_len]
                              station_δ[:, 1:old_len] .= old_station_δ[:, 1:old_len]
                            end
                            if t > tdump
                              for s = 1:numstations
                                open("$(base_name)_$(stations[s]).dat", "w") do f
                                  write(f, "t slip slip_rate shear_stress state\n")
                                  for n = 1:(nxt_ind-1)
                                    write(f, "$(station_t[n]) $(station_δ[s, n]) $(log10(abs(station_V[s, n]))) $(station_τ[s, n]) $(log10(station_θ[s, n]))\n")
                                  end
                                end
                              end
                              tdump += 100
                            end
                          end
                        end
                        Vmax
                      end, SavedValues(Float64, Float64))
  sol = solve(prob, Tsit5(); isoutofdomain=stepcheck, dt=1e3,
              atol = 1e-10, rtol = 1e-10, save_everystep=false, callback=cb)
  ψδ = sol.u[end]
  ψδ[δNp .+ (1:δNp)] = ψδ[δNp .+ (1:δNp)]
  δ = @view ψδ[δNp .+ (1:δNp)]

  for s = 1:numstations
    open("$(base_name)_$(stations[s]).dat", "w") do f
      write(f, "t slip slip_rate shear_stress state\n")
      for n = 1:(nxt_ind-1)
        write(f, "$(station_t[n]) $(station_δ[s, n]) $(log10(abs(station_V[s, n]))) $(station_τ[s, n]) $(log10(station_θ[s, n]))\n")
      end
    end
  end

  #=
  dψV = zeros(2 * δNp)
  odefun(dψV, ψδ, (), tspan[end])
  V = @view dψV[δNp .+ (1:δNp)]
  @show extrema(V)
  @show extrema(δ)

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
      Δu = u[vstarts[e]:(vstarts[e+1]-1)] # - ulinear(x,y,tspan[end])
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
  =#
  nothing
end
