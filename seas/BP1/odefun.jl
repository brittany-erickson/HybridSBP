const year_seconds = 31556926

using OrdinaryDiffEq
using DiffEqCallbacks
using Printf

function odefun(dψV, ψδ, p, t)
  RS_FAULT = 7
  VP_FAULT = 8
  reject_step = p.reject_step
  Vp = p.Vp
  lop = p.lop
  EToF = p.EToF
  EToS = p.EToS
  FToE = p.FToE
  FToLF = p.FToLF
  EToO = p.EToO
  FToB = p.FToB
  FToλstarts = p.FToλstarts
  FToδstarts = p.FToδstarts
  gδ = p.gδ
  λ = p.λ
  bλ = p.bλ
  u = p.u
  τ = p.τ
  g = p.g
  vstarts = p.vstarts
  BF = p.BF
  FbarT = p.FbarT
  locfactors = p.locfactors
  μshear = p.μshear
  RSa = p.RSa
  RSb = p.RSb
  σn = p.σn
  η = p.η
  RSV0 = p.RSV0
  τz0 = p.τz0
  RSDc = p.RSDc
  RSf0 = p.RSf0

  nelems = length(lop)
  nfaces = size(FToE, 2)
  bc_Dirichlet(lf, x, y, e, t, δ) = fill(sign(lop[e].nx[lf][1]) * t * Vp/2,
                                         size(x))
  bc_Neumann(lf, x, y, nx, ny, e, t, δ) = zeros(size(x))
  function in_jump(lf, x, y, e, t, δ)
    f = EToF[lf, e]
    if EToS[lf, e] == 1
      if EToO[lf, e]
        return -δ[FToδstarts[f]:(FToδstarts[f+1]-1)]
      else
        error("shouldn't get here")
      end
    else
      if EToO[lf, e]
        return  δ[FToδstarts[f]:(FToδstarts[f+1]-1)]
      else
        return  δ[(FToδstarts[f+1]-1):-1:FToδstarts[f]]
      end
    end
  end

  if reject_step[1]
    return
  end
  δNp = div(length(ψδ), 2)
  ψ  = @view ψδ[        (1:δNp) ]
  δ  = ψδ[ δNp .+ (1:δNp) ]

  dψ = @view dψV[       (1:δNp) ]
  V  = @view dψV[δNp .+ (1:δNp) ]

  gδ .= 0
  # Setup the interfaces/boundary data arrays
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
                FToB[EToF[:, e]], bc_Dirichlet, bc_Neumann, in_jump,
                (e, t, δ))
  end

  # Solve the global problem
  LocalToGLobalRHS!(bλ, g, gδ,  u, locfactors, FbarT, vstarts)
  λ[:] = BF \ bλ

  # Solve the local problems
  u[:] = -FbarT' * λ
  u[:] .= g .+ u
  u_calc = fill(false, nelems)
  #=
  for e = 1:nelems
    # Since we only care about he fault, we ignore the locked blocks
    F = locfactors[e]
    (x, y) = lop[e].coord
    JH = lop[e].JH

    @views u[vstarts[e]:(vstarts[e+1]-1)] = (F \
                                             u[vstarts[e]:(vstarts[e+1]-1)])
    # @views ldiv!(u[vstarts[e]:(vstarts[e+1]-1)], F,
    #              u[vstarts[e]:(vstarts[e+1]-1)])
    # @views ldiv!(F, u[vstarts[e]:(vstarts[e+1]-1)])
  end
  =#

  # Update the fault data
  τ .= 0
  for f = 1:nfaces
    e1 = FToE[1, f]
    lf1 = FToLF[1, f]
    if FToB[f] == RS_FAULT
      if !u_calc[e1]
        F = locfactors[e1]
        (x, y) = lop[e1].coord
        JH = lop[e1].JH

        @views u[vstarts[e1]:(vstarts[e1+1]-1)] =
                  (F \ u[vstarts[e1]:(vstarts[e1+1]-1)])
        u_calc[e1] = true
      end
      # Compute the shear-traction
      λrng = FToλstarts[f]:(FToλstarts[f+1]-1)
      δrng = FToδstarts[f]:(FToδstarts[f+1]-1)
      urng = vstarts[e1]:(vstarts[e1+1]-1)
      @views τ[δrng] .= μshear .* computetraction(lop[e1], lf1, u[urng],
                                                  λ[λrng], δ[δrng])
      #=
      e2 = FToE[2, f]
      lf2 = FToLF[2, f]
      τ2 = if EToO[lf2, e2]
      urng = vstarts[e2]:(vstarts[e2+1]-1)
      λrng = FToλstarts[f]:(FToλstarts[f+1]-1)
      μshear .* computetraction(lop[e2], lf2, u[urng], λ[δrng], -δ[δrng])
      else
      δrng = (FToδstarts[f+1]-1):-1:FToδstarts[f]
      λrng = (FToλstarts[f+1]-1):-1:FToλstarts[f]
      urng = vstarts[e2]:(vstarts[e2+1]-1)
      μshear .* computetraction(lop[e2], lf2, u[urng], λ[λrng], -δ[δrng])
      end
      =#

      # update velocity
      for n = 1:length(δrng)
        δn = δrng[n]
        ψn = ψ[δn]
        an = RSa[δn]
        τ[δn] = τn = τ[δn] + τz0[δn]
        if isnan(τn)
          println("τ reject")
          reject_step[1] = true
          return
        end

        VR = abs(τn / η)
        VL = -VR
        Vn = V[δn]
        obj_rs(V) = rateandstate(V, ψn, σn, τn, η, an, RSV0)
        (Vn, _, iter) = newtbndv(obj_rs, VL, VR, Vn; ftol = 1e-9,
                                 atolx = 1e-9, rtolx = 1e-9)
        if isnan(Vn) || iter < 0
          println("V reject")
          reject_step[1] = true
          return
          #error()
        end
        V[δn] = Vn

        dψ[δn] = (RSb * RSV0 / RSDc) * (exp((RSf0 - ψn) / RSb) - abs(Vn) / RSV0)
        if !isfinite(dψ[δn])
          println("ψ reject")
          dψ[δn] = 0
          reject_step[1] = true
          return
        end
      end
    elseif FToB[f] == VP_FAULT
      nx = lop[e1].nx
      V[FToδstarts[f]:(FToδstarts[f+1]-1)] .= sign(nx[lf1][1]) * Vp
      dψ[FToδstarts[f]:(FToδstarts[f+1]-1)] .= 0
    end
  end
  nothing
end


function setupfaultstations(locations, lop, FToB, FToE, FToLF, faults)
  T = eltype(locations)
  @assert size(locations, 2) == 2
  numstations = size(locations, 1)
  station_ind = zeros(Int64, numstations)
  station_face = zeros(Int64, numstations)
  nfaces = size(FToE, 2)
  for s = 1:numstations
    xs = locations[s, 1]
    ys = locations[s, 2]
    station_ind[s] = 0
    station_face[s] = 0
    d = typemax(T)
    for f = 1:nfaces
      if FToB[f] ∈ faults
        (e1, _) = FToE[:, f]
        (lf1, _) = FToLF[:, f]
        xf = lop[e1].facecoord[1][lf1]
        yf = lop[e1].facecoord[2][lf1]
        dA = hypot.(xf .- xs, yf .- ys)
        n = argmin(dA)
        if dA[n] < d
          station_ind[s] = n
          d = dA[n]
          station_face[s] = f
        end
      end
    end

    (e1, _) = FToE[:, station_face[s]]
    (lf1, _) = FToLF[:, station_face[s]]
    xf = lop[e1].facecoord[1][lf1][station_ind[s]]
    yf = lop[e1].facecoord[2][lf1][station_ind[s]]
  end
  return (ind=station_ind, face=station_face,
          xs = locations[:, 1],
          ys = locations[:, 2],
          tnext=zeros(T, 1),
          tdump=zeros(T, 1),
          t=Array{T, 1}(undef, 0),
          data=ntuple(n->(V=Array{T, 1}(undef, 0),
                          τ=Array{T, 1}(undef, 0),
                          θ=Array{T, 1}(undef, 0),
                          δ=Array{T, 1}(undef, 0)),
                      numstations))
end

function savefaultstation(ψδ, t, i, stations, FToδstarts, p, base_name="",
                         tdump=100)
  Vmax = 0.0
  if isdefined(i, :fsallast)
    δNp = div(length(ψδ), 2)
    dψV = i.fsallast
    V  = @view dψV[δNp .+ (1:δNp) ]
    Vmax = maximum(abs.(extrema(V)))
    tlast = length(stations.t) > 0 ? stations.t[end] : -2year_seconds
    tnext = tlast + (Vmax > 1e-3 ? 0.1 : year_seconds)
    if (t >= tnext)
      ψ  = @view ψδ[        (1:δNp) ]
      δ  = @view ψδ[ δNp .+ (1:δNp) ]
      dψ = @view dψV[       (1:δNp) ]
      tlast = tnext
      @show (t/year_seconds, Vmax)
      push!(stations.t, t)
      numstations = length(stations.data)
      for s = 1:numstations
        f = stations.face[s]
        n = stations.ind[s] + FToδstarts[f] - 1
        push!(stations.data[s].V, V[n])
        push!(stations.data[s].θ, (p.RSDc * exp((ψ[n] - p.RSf0) / p.RSb) /
                                   p.RSV0))
        push!(stations.data[s].δ, δ[n])
        push!(stations.data[s].τ, p.τ[n] - p.η * V[n])
      end
      if length(stations.t) == 1 ||
        ceil((stations.t[end] / tdump)) > ceil((stations.t[end-1] / tdump))
        println("saving data with basename = $base_name")
        for s = 1:numstations
          open("$(base_name)$(stations.xs[s])_$(stations.ys[s]).dat", "w") do f
            write(f, "t slip slip_rate shear_stress state\n")
            t = stations.t
            δ = stations.data[s].δ
            V = stations.data[s].V
            θ = stations.data[s].θ
            τ = stations.data[s].τ
            for n = 1:length(t)
              write(f, "$(t[n]) $(δ[n]) $(log10(abs(V[n]))) $(τ[n]) $(log10(θ[n]))\n")
            end
          end
        end
      end
    end
  end
  Vmax
end
