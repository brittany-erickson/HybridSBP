if VERSION == v"0.6.3"
  eigen = eig
  macro plotting(ex)
    return :($(esc(ex)))
  end
else
  macro plotting(ex)
  end
end

@plotting let
  using Plots
  pyplot()
end

include("diagonal_sbp.jl")

using Compat
import Compat: range, undef
using Compat.SparseArrays

⊗ = (A,B) -> kron(A, B)

function locoperator(p, Nr, Ns)
  Nrp = Nr + 1
  Nsp = Ns + 1
  Np = Nrp * Nsp

  (D2r, Sr0, SrN, HIr, Hr, r) = diagonal_sbp_D2(p, Nr)
  (D2s, Ss0, SsN, HIs, Hs, s) = diagonal_sbp_D2(p, Ns)

  Ir = sparse(1.0I, Nrp, Nrp)
  Is = sparse(1.0I, Nsp, Nsp)

  (r, s) = (ones(Nsp) ⊗ r, s ⊗ ones(Nrp))

  Ar = SrN - Sr0- Hr * D2r
  As = SsN - Ss0- Hs * D2s

  E0r = sparse([     1], [     1], [1], Nr + 1, Nr + 1)
  ENr = sparse([Nr + 1], [Nr + 1], [1], Nr + 1, Nr + 1)
  E0s = sparse([     1], [     1], [1], Ns + 1, Ns + 1)
  ENs = sparse([Ns + 1], [Ns + 1], [1], Ns + 1, Ns + 1)

  e0r = sparse([     1], [1], [1], Nr + 1, 1)
  eNr = sparse([Nr + 1], [1], [1], Nr + 1, 1)
  e0s = sparse([     1], [1], [1], Ns + 1, 1)
  eNs = sparse([Ns + 1], [1], [1], Ns + 1, 1)

  hr = 2 / Nr
  hs = 2 / Ns

  τr = 2 * max(Nr, Ns)
  τs = 2 * max(Nr, Ns)
  M = (Hs ⊗ (Ar - SrN + Sr0 - SrN' + Sr0' + τr * (E0r + ENr))) +
      ((As - SsN + Ss0 - SsN' + Ss0' + τs * (E0s + ENs)) ⊗ Hr)

  er0 = sparse([1  ], [1], [1], Nrp, 1)
  erN = sparse([Nrp], [1], [1], Nrp, 1)
  es0 = sparse([1  ], [1], [1], Nsp, 1)
  esN = sparse([Nsp], [1], [1], Nsp, 1)

  L1 = (Is ⊗ er0')
  L2 = (Is ⊗ erN')
  L3 = (es0' ⊗ Ir)
  L4 = (esN' ⊗ Ir)

  F1 =  (Hs ⊗ (Sr0' * er0)) + ((τs * Hs) ⊗ er0)
  F2 = -(Hs ⊗ (SrN' * erN)) + ((τs * Hs) ⊗ erN)
  F3 =  ((Ss0' * es0) ⊗ Hr) + (es0 ⊗ (τr * Hr))
  F4 = -((SsN' * esN) ⊗ Hr) + (esN ⊗ (τr * Hr))

  (M, (F1, F2, F3, F4), (L1, L2, L3, L4), (r, s), Hs ⊗ Hr)
end

function transfinite_blend(α1, α2, α3, α4, r, s)
  # +---4---+
  # |       |
  # 1       2
  # |       |
  # +---3---+
  @assert α1(-1) ≈ α3(-1)
  @assert α2(-1) ≈ α3( 1)
  @assert α1( 1) ≈ α4(-1)
  @assert α2( 1) ≈ α4( 1)

  x = (1 .+ r) .* α2(s)/2 + (1 .- r) .* α1(s)/2 +
      (1 .+ s) .* α4(r)/2 + (1 .- s) .* α3(r)/2 -
     ((1 .+ r) .* (1 .+ s) .* α2( 1) +
      (1 .- r) .* (1 .+ s) .* α1( 1) +
      (1 .+ r) .* (1 .- s) .* α2(-1) +
      (1 .- r) .* (1 .- s) .* α1(-1)) / 4
end


function locoperator(p, Nr, Ns, xf, yf; pm = p+2)
  Nrp = Nr + 1
  Nsp = Ns + 1
  Np = Nrp * Nsp

  (DrM, ~, ~, ~) = diagonal_sbp_D1(pm, Nr; xc = (-1,1))
  (DsM, ~, ~, ~) = diagonal_sbp_D1(pm, Ns; xc = (-1,1))

  (Dr, HrI, Hr, r) = diagonal_sbp_D1(p, Nr; xc = (-1,1))
  (Ds, HsI, Hs, s) = diagonal_sbp_D1(p, Ns; xc = (-1,1))

  Ir = sparse(1.0I, Nrp, Nrp)
  Is = sparse(1.0I, Nsp, Nsp)

  Qr = Hr * Dr
  Qs = Hs * Ds

  (r, s) = (ones(Nsp) ⊗ r, s ⊗ ones(Nrp))
  (x, y) = (xf(r, s), yf(r, s))

  # Compute the metric terms
  xr = (Is ⊗ DrM) * x
  xs = (DsM ⊗ Ir) * x
  yr = (Is ⊗ DrM) * y
  ys = (DsM ⊗ Ir) * y

  J = xr .* ys - xs .* yr
  @assert minimum(J) > 0

  rx =  ys ./ J
  sx = -yr ./ J
  ry = -xs ./ J
  sy =  xr ./ J

  crr = J .* (rx .* rx + ry .* ry)
  crs = csr = J .* (sx .* rx + sy .* ry)
  css = J .* (sx .* sx + sy .* sy)

  IArr = Array{Int64,1}(undef,0)
  JArr = Array{Int64,1}(undef,0)
  VArr = Array{Float64,1}(undef,0)
  ISr0 = Array{Int64,1}(undef,0)
  JSr0 = Array{Int64,1}(undef,0)
  VSr0 = Array{Float64,1}(undef,0)
  ISrN = Array{Int64,1}(undef,0)
  JSrN = Array{Int64,1}(undef,0)
  VSrN = Array{Float64,1}(undef,0)
  for j = 1:Nsp
    rng = (j-1) * Nrp .+ (1:Nrp)
    (~, S0e, SNe, ~, ~, Ae, ~) = variable_diagonal_sbp_D2(p, Nr, crr[rng])
    (Ie, Je, Ve) = findnz(Ae)
    IArr = [IArr;Ie .+ (j-1) * Nrp]
    JArr = [JArr;Je .+ (j-1) * Nrp]
    VArr = [VArr;Hs[j,j] * Ve]

    (Ie, Je, Ve) = findnz(S0e)
    ISr0 = [ISr0;Ie .+ (j-1) * Nrp]
    JSr0 = [JSr0;Je .+ (j-1) * Nrp]
    VSr0 = [VSr0; Hs[j,j] * Ve]

    (Ie, Je, Ve) = findnz(SNe)
    ISrN = [ISrN;Ie .+ (j-1) * Nrp]
    JSrN = [JSrN;Je .+ (j-1) * Nrp]
    VSrN = [VSrN; Hs[j,j] * Ve]
  end
  Arr = sparse(IArr, JArr, VArr, Np, Np)
  Sr0 = sparse(ISr0, JSr0, VSr0, Np, Np)
  SrN = sparse(ISrN, JSrN, VSrN, Np, Np)
  @assert Arr ≈ Arr'
  (D2, S0, SN, ~, ~, ~) = diagonal_sbp_D2(p, Nr)
  #= affine mesh test
  Ar = SN - S0 - Hr * D2
  @assert Arr ≈ Hs ⊗ Ar
  =#
  @assert Sr0 ≈ ((sparse(Diagonal(crr[1   .+ Nrp*(0:Ns)])) * Hs) ⊗ S0)
  @assert SrN ≈ ((sparse(Diagonal(crr[Nrp .+ Nrp*(0:Ns)])) * Hs) ⊗ SN)

  IAss = Array{Int64,1}(undef,0)
  JAss = Array{Int64,1}(undef,0)
  VAss = Array{Float64,1}(undef,0)
  ISs0 = Array{Int64,1}(undef,0)
  JSs0 = Array{Int64,1}(undef,0)
  VSs0 = Array{Float64,1}(undef,0)
  ISsN = Array{Int64,1}(undef,0)
  JSsN = Array{Int64,1}(undef,0)
  VSsN = Array{Float64,1}(undef,0)
  for i = 1:Nrp
    rng = i .+ Nrp * (0:Ns)
    (~, S0e, SNe, ~, ~, Ae, ~) = variable_diagonal_sbp_D2(p, Ns, css[rng])

    (Ie, Je, Ve) = findnz(Ae)
    IAss = [IAss;i .+ Nrp * (Ie .- 1)]
    JAss = [JAss;i .+ Nrp * (Je .- 1)]
    VAss = [VAss;Hr[i,i] * Ve]

    (Ie, Je, Ve) = findnz(S0e)
    ISs0 = [ISs0;i .+ Nrp * (Ie .- 1)]
    JSs0 = [JSs0;i .+ Nrp * (Je .- 1)]
    VSs0 = [VSs0;Hr[i,i] * Ve]

    (Ie, Je, Ve) = findnz(SNe)
    ISsN = [ISsN;i .+ Nrp * (Ie .- 1)]
    JSsN = [JSsN;i .+ Nrp * (Je .- 1)]
    VSsN = [VSsN;Hr[i,i] * Ve]
  end
  Ass = sparse(IAss, JAss, VAss, Np, Np)
  Ss0 = sparse(ISs0, JSs0, VSs0, Np, Np)
  SsN = sparse(ISsN, JSsN, VSsN, Np, Np)
  @assert Ass ≈ Ass'
  (D2, S0, SN, ~, ~, ~) = diagonal_sbp_D2(p, Ns)
  #= affine mesh test
  As = SN - S0 - Hs * D2
  @assert Ass ≈ As ⊗ Hr
  =#
  @assert Ss0 ≈ (S0 ⊗ (Hr * sparse(Diagonal(css[1:Nrp]))))
  @assert SsN ≈ (SN ⊗ (Hr * sparse(Diagonal(css[Nrp*Ns .+ (1:Nrp)]))))

  Asr = (sparse(Qs') ⊗ Ir) * sparse(Diagonal(crs)) * (Is ⊗ Qr)
  Ars = (Is ⊗ sparse(Qr')) * sparse(Diagonal(csr)) * (Qs ⊗ Ir)

  A = Arr + Ass + Ars + Asr

  Er0 = sparse([1], [1], [1], Nrp, Nrp)
  ErN = sparse([Nrp], [Nrp], [1], Nrp, Nrp)
  Es0 = sparse([1], [1], [1], Nsp, Nsp)
  EsN = sparse([Nsp], [Nsp], [1], Nsp, Nsp)

  crs0 = sparse(Diagonal(crs[1:Nrp]))
  crsN = sparse(Diagonal(crs[Nrp*Ns .+ (1:Nrp)]))
  csr0 = sparse(Diagonal(csr[1   .+ Nrp*(0:Ns)]))
  csrN = sparse(Diagonal(csr[Nrp .+ Nrp*(0:Ns)]))

  τ1 = sparse(Diagonal(2*max(Nr, Ns)*ones(Nsp)))
  τ2 = sparse(Diagonal(2*max(Nr, Ns)*ones(Nsp)))
  τ3 = sparse(Diagonal(2*max(Nr, Ns)*ones(Nrp)))
  τ4 = sparse(Diagonal(2*max(Nr, Ns)*ones(Nrp)))

  er0 = sparse([1  ], [1], [1], Nrp, 1)
  erN = sparse([Nrp], [1], [1], Nrp, 1)
  es0 = sparse([1  ], [1], [1], Nsp, 1)
  esN = sparse([Nsp], [1], [1], Nsp, 1)

  F1 =  Sr0' * (Is ⊗ er0) + ((Qs' * csr0) ⊗ er0) + ((τ1 * Hs) ⊗ er0)
  F2 = -SrN' * (Is ⊗ erN) - ((Qs' * csrN) ⊗ erN) + ((τ2 * Hs) ⊗ erN)
  F3 =  Ss0' * (es0 ⊗ Ir) + (es0 ⊗ Qr' * crs0) + (es0 ⊗ (τ3 * Hr))
  F4 = -SsN' * (esN ⊗ Ir) - (esN ⊗ Qr' * crsN) + (esN ⊗ (τ4 * Hr))

  # TODO: Check signs on Q terms (and update write up with correct signs)
  B1 =  (Sr0 + Sr0') + ((csr0 * Qs + Qs' * csr0) ⊗ Er0) + ((τ1 * Hs) ⊗ Er0)
  B2 = -(SrN + SrN') - ((csrN * Qs + Qs' * csrN) ⊗ ErN) + ((τ2 * Hs) ⊗ ErN)
  B3 =  (Ss0 + Ss0') + (Es0 ⊗ (crs0 * Qr + Qr' * crs0)) + (Es0 ⊗ (τ3 * Hr))
  B4 = -(SsN + SsN') - (EsN ⊗ (crsN * Qr + Qr' * crsN)) + (EsN ⊗ (τ4 * Hr))

  @assert B1 ≈ F1 * (Is ⊗ er0') + (Is ⊗ er0) * F1' - ((τ1 * Hs) ⊗ Er0)
  @assert B2 ≈ F2 * (Is ⊗ erN') + (Is ⊗ erN) * F2' - ((τ2 * Hs) ⊗ ErN)
  @assert B3 ≈ F3 * (es0' ⊗ Ir) + (es0 ⊗ Ir) * F3' - (Es0 ⊗ (τ3 * Hr))
  @assert B4 ≈ F4 * (esN' ⊗ Ir) + (esN ⊗ Ir) * F4' - (EsN ⊗ (τ4 * Hr))

  M = A + B1 + B2 + B3 + B4

  L1 = (Is ⊗ er0')
  L2 = (Is ⊗ erN')
  L3 = (es0' ⊗ Ir)
  L4 = (esN' ⊗ Ir)
  (x1, y1) = (L1 * x, L1 * y)
  (x2, y2) = (L2 * x, L2 * y)
  (x3, y3) = (L3 * x, L3 * y)
  (x4, y4) = (L4 * x, L4 * y)

  #=
  @plotting let
    plot(x1, y1)
    plot!(x2, y2)
    plot!(x3, y3)
    plot!(x4, y4)
    display(plot!())
  end
  =#

  # (E, V) = eigen(Matrix(M))
  # println((minimum(E), maximum(E)))
  (M, (F1, F2, F3, F4), (L1, L2, L3, L4), (x, y), Diagonal(J) * (Hs ⊗ Hr))
end

let
  p = 4
  Nr0 = 20
  Ns0 = 21

  x1 = (s)-> sin.(π*s)/8 .- 1
  y1 = (s)-> s

  x2 = (s)->sin.(2*π*s)/8 .+ 1
  y2 = (s)-> s

  x3 = (r)-> r
  y3 = (r)-> sin.(2*π*r)/8 .- 1

  x4 = (r)->r
  y4 = (r)->-sin.(π*r)/8 .+ 1

  xt = (r,s)->transfinite_blend(x1, x2, x3, x4, r, s)
  yt = (r,s)->transfinite_blend(y1, y2, y3, y4, r, s)

  # xt = (r,s)->r
  # yt = (r,s)->s

  @plotting let
    r = Compat.range(-1, stop=1, length=101)
    s = ones(101)

    plot(xt(-s, r), yt(-s, r))
    plot!(xt(s, r), yt(s, r))
    plot!(xt(r, -s), yt(r, -s))
    plot!(xt(r, s), yt(r, s))

    r = Compat.range(-1, stop=1, length=20) * ones(1, 20)
    s = (Compat.range(-1, stop=1, length=20) * ones(1, 20))'
    plot!(xt(r, s), yt(r, s))
    plot!(xt(r', s'), yt(r', s'))
    display(plot!())
  end

  ϵ = zeros(5)
  (kx, ky) = (π, π)
  vex = (x,y) -> cos.(kx * x) .* cosh.(ky * y)
  for lvl = 1:length(ϵ)
    Nr = Nr0 * (2^(lvl-1))
    Ns = Ns0 * (2^(lvl-1))
    (M, (F1, F2, F3, F4), (L1, L2, L3, L4), (x, y), H) =
      locoperator(p, Nr, Ns, xt, yt)
    #=
    (M, (F1, F2, F3, F4), (L1, L2, L3, L4), (x, y), H) =
      locoperator(p, Nr, Ns)
    =#
    v = vex(x,y)
    v1 = L1 * v
    v2 = L2 * v
    v3 = L3 * v
    v4 = L4 * v
    g = F1 * v1 + F2 * v2 + F3 * v3 + F4 * v4
    u = M \ g
    Δ = u - v
    ϵ[lvl] = sqrt(Δ' * H * Δ)
    println("level = ", lvl, " :: error = ", ϵ[lvl])
    @plotting let
      if lvl == 1
        l = @layout grid(2,1)
        plot(reshape(x, Nr+1, Ns+1),
             reshape(y, Nr+1, Ns+1),
        [reshape(Δ, Nr+1, Ns+1), reshape(v, Nr+1, Ns+1)],
        st = [:surface, :surface], layout=l)
        display(plot!())
      end
    end
  end
  println(ϵ)
  println((log.(ϵ[1:end-1]) - log.(ϵ[2:end])) / log(2))
  nothing
end
