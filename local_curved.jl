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


function locoperator(p, Nr, Ns, xf, yf)
  Nrp = Nr + 1
  Nsp = Ns + 1
  Np = Nrp * Nsp

  (Dr, HrI, Hr, r) = diagonal_sbp_D1(p, Nr; xc = (-1,1))
  (Ds, HsI, Hs, s) = diagonal_sbp_D1(p, Ns; xc = (-1,1))

  Ir = sparse(1.0I, Nrp, Nrp)
  Is = sparse(1.0I, Nsp, Nsp)

  Qr = Hr * Dr
  Qs = Hs * Ds

  (r, s) = (ones(Nsp) ⊗ r, s ⊗ ones(Nrp))
  (x, y) = (xf(r, s), yf(r, s))

  # Compute the metric terms
  xr = (Is ⊗ Dr) * x
  xs = (Ds ⊗ Ir) * x
  yr = (Is ⊗ Dr) * y
  ys = (Ds ⊗ Ir) * y

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
    VArr = [VArr;Ve]

    (Ie, Je, Ve) = findnz(S0e)
    ISr0 = [ISr0;Ie .+ (j-1) * Nrp]
    JSr0 = [JSr0;Je .+ (j-1) * Nrp]
    VSr0 = [VSr0;Ve]

    (Ie, Je, Ve) = findnz(SNe)
    ISrN = [ISrN;Ie .+ (j-1) * Nrp]
    JSrN = [JSrN;Je .+ (j-1) * Nrp]
    VSrN = [VSrN;Ve]
  end
  Arr = sparse(IArr, JArr, VArr, Np, Np)
  Sr0 = sparse(ISr0, JSr0, VSr0, Np, Np)
  SrN = sparse(ISrN, JSrN, VSrN, Np, Np)
  @assert Arr ≈ Arr'
  (~, S0, SN, ~, ~, ~) = diagonal_sbp_D2(p, Nr)
  @assert Sr0 ≈ (sparse(Diagonal(crr[1   .+ Nrp*(0:Ns)])) ⊗ S0)
  @assert SrN ≈ (sparse(Diagonal(crr[Nrp .+ Nrp*(0:Ns)])) ⊗ SN)

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
    VAss = [VAss;Ve]

    (Ie, Je, Ve) = findnz(S0e)
    ISs0 = [ISs0;i .+ Nrp * (Ie .- 1)]
    JSs0 = [JSs0;i .+ Nrp * (Je .- 1)]
    VSs0 = [VSs0;Ve]

    (Ie, Je, Ve) = findnz(SNe)
    ISsN = [ISsN;i .+ Nrp * (Ie .- 1)]
    JSsN = [JSsN;i .+ Nrp * (Je .- 1)]
    VSsN = [VSsN;Ve]
  end
  Ass = sparse(IAss, JAss, VAss, Np, Np)
  Ss0 = sparse(ISs0, JSs0, VSs0, Np, Np)
  SsN = sparse(ISsN, JSsN, VSsN, Np, Np)
  @assert Ass ≈ Ass'
  (~, S0, SN, ~, ~, r1d) = diagonal_sbp_D2(p, Ns)
  @assert Ss0 ≈ (S0 ⊗ sparse(Diagonal(css[1:Nrp])))
  @assert SsN ≈ (SN ⊗ sparse(Diagonal(css[Nrp*Ns .+ (1:Nrp)])))

  Ars = (Qs' ⊗ Ir) * sparse(Diagonal(crs)) * (Is ⊗ Qr)
  Asr = (Is ⊗ Qr') * sparse(Diagonal(csr)) * (Qs ⊗ Ir)

  A = Arr + Ass + Ars + Asr

  Er0 = sparse([1], [1], [1], Nrp, Nrp)
  ErN = sparse([Nrp], [Nrp], [1], Nrp, Nrp)
  Es0 = sparse([1], [1], [1], Nsp, Nsp)
  EsN = sparse([Nsp], [Nsp], [1], Nsp, Nsp)

  crs0 = sparse(Diagonal(crs[1:Nrp]))
  crsN = sparse(Diagonal(crs[Nrp*Ns .+ (1:Nrp)]))
  csr0 = sparse(Diagonal(csr[1   .+ Nrp*(0:Ns)]))
  csrN = sparse(Diagonal(csr[Nrp .+ Nrp*(0:Ns)]))

  τ1 = sparse(Diagonal(10*ones(Nrp)))
  τ2 = sparse(Diagonal(10*ones(Nrp)))
  τ3 = sparse(Diagonal(10*ones(Nsp)))
  τ4 = sparse(Diagonal(10*ones(Nsp)))

  # TODO: Check signs on Q terms (and update write up with correct signs)
  B1 =  (Ss0 + Ss0') * (Is ⊗ Hr) + (Es0 ⊗ (crs0 * Qr + Qr' * crs0)) + (Es0 ⊗ (τ1 * Hr))
  B2 = -(SsN + SsN') * (Is ⊗ Hr) - (EsN ⊗ (crsN * Qr + Qr' * crsN)) + (EsN ⊗ (τ2 * Hr))
  B3 =  (Sr0 + Sr0') * (Hs ⊗ Ir) + ((csr0 * Qs + Qs' * csr0) ⊗ Er0) + ((τ3 * Hs) ⊗ Er0)
  B4 = -(SrN + SrN') * (Hs ⊗ Ir) - ((csrN * Qs + Qs' * csrN) ⊗ ErN) + ((τ4 * Hs) ⊗ ErN)

  M = A + B1 + B2 + B3 + B4

  er0 = sparse([1  ], [1], [1], Nrp, 1)
  erN = sparse([Nrp], [1], [1], Nrp, 1)
  es0 = sparse([1  ], [1], [1], Nsp, 1)
  esN = sparse([Nsp], [1], [1], Nsp, 1)

  F1 =  Ss0' * (es0 ⊗ Hr) + (es0 ⊗ Qr' * crs0) + (es0 ⊗ τ1 * Hr)
  F2 = -SsN' * (esN ⊗ Hr) - (esN ⊗ Qr' * crsN) + (esN ⊗ τ2 * Hr)
  F3 =  (Hs ⊗ er0) * Sr0' #+ (Qs' * csr0) ⊗ er0) + (τ3 * Hs ⊗ er0)
  F4 = -(Hs ⊗ erN) * SrN' #+ (Qs' * csrN) ⊗ erN) + (τ4 * Hs ⊗ erN)
  println(size(F1))
  println(size(F2))
  println(size(F3))
  println(size(F4))

  (E, V) = eigen(Matrix(M))
  println((minimum(E), maximum(E)))
end

let
  p = 4
  Nr = 20
  Ns = 30

  x1 = (s)-> sin.(π*s)/8 .- 1
  y1 = (s)-> s

  x2 = (s)->sin.(2*π*s)/8 .+ 1
  y2 = (s)-> s

  x3 = (r)-> r
  y3 = (r)-> sin.(2*π*r)/8 .- 1

  x4 = (r)->r
  y4 = (r)->-sin.(π*r)/8 .+ 1

  x = (r,s)->transfinite_blend(x1, x2, x3, x4, r, s)
  y = (r,s)->transfinite_blend(y1, y2, y3, y4, r, s)

  @plotting let
    r = Compat.range(-1, stop=1, length=101)
    s = ones(101)

    plot(x(-s, r), y(-s, r))
    plot!(x(s, r), y(s, r))
    plot!(x(r, -s), y(r, -s))
    plot!(x(r, s), y(r, s))

    r = Compat.range(-1, stop=1, length=20) * ones(1, 20)
    s = (Compat.range(-1, stop=1, length=20) * ones(1, 20))'
    plot!(x(r, s), y(r, s))
    plot!(x(r', s'), y(r', s'))
    display(plot!())
  end

  locoperator(p, Nr, Ns, x, y)
  nothing
end
