const year_seconds = 31556926

using OrdinaryDiffEq
using DiffEqCallbacks
using Printf
using Plots

function odefun(dψV, ψδ, p, t)
  reject_step = p.reject_step
  Vp = p.Vp
  lop = p.lop
  F = p.F
  u = p.u
  Δτ = p.Δτ
  ge = p.ge
  μshear = p.μshear
  RSa = p.RSa
  RSb = p.RSb
  σn = p.σn
  η = p.η
  RSV0 = p.RSV0
  τz0 = p.τz0
  RSDc = p.RSDc
  RSf0 = p.RSf0
  LFtoB = p.LFtoB


  if reject_step[1]
    return
  end
  δNp = div(length(ψδ), 2)
  ψ  = @view ψδ[        (1:δNp) ]
  δ  = ψδ[ δNp .+ (1:δNp) ]


  bc_Dirichlet = (lf, x, y, e) -> (2-lf)*(δ ./ 2) + (lf-1)*fill(t * Vp/2, size(x))
  bc_Neumann   = (lf, x, y, nx, ny, e) -> zeros(size(x))


  # solve for displacements everywhere in domain
  e = 1
  locbcarray_mod!(ge, lop[e], LFtoB, bc_Dirichlet, bc_Neumann,(e))
  u[:] = F \ ge


  # set up rates of change for  state and slip
  dψ = @view dψV[       (1:δNp) ]
  V  = @view dψV[δNp .+ (1:δNp) ]

  dψ .= 0
  V  .= 0


  # Update the fault data
  Δτ .= 0
  lf1 = 1


  Δτ .= -μshear .* computetraction_mod(lop[1], lf1, u, δ)

  #δex = Vp*t/2
  #τ_ex = μshear * (Vp * t/2 -δex/2)/20
  Vex = zeros(δNp)
  #yf = lop[1].facecoord[2][1]
  #Vtest = Vp/2
  # update velocity
#@show t./31556926

  for n = 1:δNp
    ψn = ψ[n]
    an = RSa[n]

    τn = Δτ[n] + τz0
    if isnan(τn)
      println("τ reject")
      reject_step[1] = true
      return
    end

    VR = abs(τn / η)
    VL = -VR
    Vn = V[n]
    obj_rs(V) = rateandstate(V, ψn, σn, τn, η, an, RSV0)
    (Vn, _, iter) = newtbndv(obj_rs, VL, VR, Vn; ftol = 1e-9,
                                 atolx = 1e-9, rtolx = 1e-9)

    Vex[n] = 2*RSV0*sinh(τn/(σn * an)) * exp(-ψn / an)



    if isnan(Vn) || iter < 0
      println("V reject")
      reject_step[1] = true
      return
          #error()
    end
    V[n] = Vn



    dψ[n] = (RSb * RSV0 / RSDc) * (exp((RSf0 - ψn) / RSb) - abs(Vn) / RSV0)
    if !isfinite(dψ[n])
      println("ψ reject")
      dψ[n] = 0
      reject_step[1] = true
      return
    end
  end
  #@show δ[end]
  #sleep(2)
  #@show [maximum(V) minimum(ψ) maximum(τ)]
  #pyplot()

  #display(plot(yf, δ))
  #sleep(1)
 #@show [minimum(V) maximum(V) t/31556926]
 #sleep(.1)
#  @show t/31556926
#sleep(1)
  nothing
end
