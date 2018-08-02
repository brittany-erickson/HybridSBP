using Plots
pyplot()

if VERSION <= v"0.6.999999"
  eigen = eig
  macro isdefined(s::Symbol)
    return isdefined(s)
  end
end
include("diagonal_sbp.jl")
if !@isdefined ⊗
  const ⊗ = (A,B) -> kron(A, B)
end

#{{{ locoperator
function locoperator(p, Nr, Ns, xf, yf; pm = p+2, LFToB = [])
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
  QrT = sparse(transpose(Qr))
  Qs = Hs * Ds
  QsT = sparse(transpose(Qs))

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

  ISr0 = Array{Int64,1}(undef,0)
  JSr0 = Array{Int64,1}(undef,0)
  VSr0 = Array{Float64,1}(undef,0)
  ISrN = Array{Int64,1}(undef,0)
  JSrN = Array{Int64,1}(undef,0)
  VSrN = Array{Float64,1}(undef,0)

  (~, S0e, SNe, ~, ~, Ae, ~) = variable_diagonal_sbp_D2(p, Nr, rand(Nrp))
  IArr = Array{Int64,1}(undef,Nsp * length(Ae.nzval))
  JArr = Array{Int64,1}(undef,Nsp * length(Ae.nzval))
  VArr = Array{Float64,1}(undef,Nsp * length(Ae.nzval))
  stArr = 0

  ISr0 = Array{Int64,1}(undef,Nsp * length(S0e.nzval))
  JSr0 = Array{Int64,1}(undef,Nsp * length(S0e.nzval))
  VSr0 = Array{Float64,1}(undef,Nsp * length(S0e.nzval))
  stSr0 = 0

  ISrN = Array{Int64,1}(undef,Nsp * length(SNe.nzval))
  JSrN = Array{Int64,1}(undef,Nsp * length(SNe.nzval))
  VSrN = Array{Float64,1}(undef,Nsp * length(SNe.nzval))
  stSrN = 0
  for j = 1:Nsp
    rng = (j-1) * Nrp .+ (1:Nrp)
    (~, S0e, SNe, ~, ~, Ae, ~) = variable_diagonal_sbp_D2(p, Nr, crr[rng])
    (Ie, Je, Ve) = findnz(Ae)
    IArr[stArr .+ (1:length(Ve))] = Ie .+ (j-1) * Nrp
    JArr[stArr .+ (1:length(Ve))] = Je .+ (j-1) * Nrp
    VArr[stArr .+ (1:length(Ve))] = Hs[j,j] * Ve
    stArr += length(Ve)

    (Ie, Je, Ve) = findnz(S0e)
    ISr0[stSr0 .+ (1:length(Ve))] = Ie .+ (j-1) * Nrp
    JSr0[stSr0 .+ (1:length(Ve))] = Je .+ (j-1) * Nrp
    VSr0[stSr0 .+ (1:length(Ve))] =  Hs[j,j] * Ve
    stSr0 += length(Ve)

    (Ie, Je, Ve) = findnz(SNe)
    ISrN[stSrN .+ (1:length(Ve))] = Ie .+ (j-1) * Nrp
    JSrN[stSrN .+ (1:length(Ve))] = Je .+ (j-1) * Nrp
    VSrN[stSrN .+ (1:length(Ve))] =  Hs[j,j] * Ve
    stSrN += length(Ve)
  end
  Arr = sparse(IArr[1:stArr], JArr[1:stArr], VArr[1:stArr], Np, Np)
  Sr0 = sparse(ISr0[1:stSr0], JSr0[1:stSr0], VSr0[1:stSr0], Np, Np)
  SrN = sparse(ISrN[1:stSrN], JSrN[1:stSrN], VSrN[1:stSrN], Np, Np)
  Sr0T = sparse(JSr0[1:stSr0], ISr0[1:stSr0], VSr0[1:stSr0], Np, Np)
  SrNT = sparse(JSrN[1:stSrN], ISrN[1:stSrN], VSrN[1:stSrN], Np, Np)
  # @assert Arr ≈ Arr'
  (D2, S0, SN, ~, ~, ~) = diagonal_sbp_D2(p, Nr)
  #= affine mesh test
  Ar = SN - S0 - Hr * D2
  @assert Arr ≈ Hs ⊗ Ar
  =#
  # @assert Sr0 ≈ ((sparse(Diagonal(crr[1   .+ Nrp*(0:Ns)])) * Hs) ⊗ S0)
  # @assert SrN ≈ ((sparse(Diagonal(crr[Nrp .+ Nrp*(0:Ns)])) * Hs) ⊗ SN)

  (~, S0e, SNe, ~, ~, Ae, ~) = variable_diagonal_sbp_D2(p, Ns, rand(Nsp))
  IAss = Array{Int64,1}(undef,Nrp * length(Ae.nzval))
  JAss = Array{Int64,1}(undef,Nrp * length(Ae.nzval))
  VAss = Array{Float64,1}(undef,Nrp * length(Ae.nzval))
  stAss = 0

  ISs0 = Array{Int64,1}(undef,Nrp * length(S0e.nzval))
  JSs0 = Array{Int64,1}(undef,Nrp * length(S0e.nzval))
  VSs0 = Array{Float64,1}(undef,Nrp * length(S0e.nzval))
  stSs0 = 0

  ISsN = Array{Int64,1}(undef,Nrp * length(SNe.nzval))
  JSsN = Array{Int64,1}(undef,Nrp * length(SNe.nzval))
  VSsN = Array{Float64,1}(undef,Nrp * length(SNe.nzval))
  stSsN = 0
  for i = 1:Nrp
    rng = i .+ Nrp * (0:Ns)
    (~, S0e, SNe, ~, ~, Ae, ~) = variable_diagonal_sbp_D2(p, Ns, css[rng])

    (Ie, Je, Ve) = findnz(Ae)
    IAss[stAss .+ (1:length(Ve))] = i .+ Nrp * (Ie .- 1)
    JAss[stAss .+ (1:length(Ve))] = i .+ Nrp * (Je .- 1)
    VAss[stAss .+ (1:length(Ve))] = Hr[i,i] * Ve
    stAss += length(Ve)

    (Ie, Je, Ve) = findnz(S0e)
    ISs0[stSs0 .+ (1:length(Ve))] = i .+ Nrp * (Ie .- 1)
    JSs0[stSs0 .+ (1:length(Ve))] = i .+ Nrp * (Je .- 1)
    VSs0[stSs0 .+ (1:length(Ve))] = Hr[i,i] * Ve
    stSs0 += length(Ve)

    (Ie, Je, Ve) = findnz(SNe)
    ISsN[stSsN .+ (1:length(Ve))] = i .+ Nrp * (Ie .- 1)
    JSsN[stSsN .+ (1:length(Ve))] = i .+ Nrp * (Je .- 1)
    VSsN[stSsN .+ (1:length(Ve))] = Hr[i,i] * Ve
    stSsN += length(Ve)
  end
  Ass = sparse(IAss[1:stAss], JAss[1:stAss], VAss[1:stAss], Np, Np)
  Ss0 = sparse(ISs0[1:stSs0], JSs0[1:stSs0], VSs0[1:stSs0], Np, Np)
  SsN = sparse(ISsN[1:stSsN], JSsN[1:stSsN], VSsN[1:stSsN], Np, Np)
  Ss0T = sparse(JSs0[1:stSs0], ISs0[1:stSs0], VSs0[1:stSs0], Np, Np)
  SsNT = sparse(JSsN[1:stSsN], ISsN[1:stSsN], VSsN[1:stSsN], Np, Np)
  # @assert Ass ≈ Ass'
  (D2, S0, SN, ~, ~, ~) = diagonal_sbp_D2(p, Ns)
  #= affine mesh test
  As = SN - S0 - Hs * D2
  @assert Ass ≈ As ⊗ Hr
  =#
  # @assert Ss0 ≈ (S0 ⊗ (Hr * sparse(Diagonal(css[1:Nrp]))))
  # @assert SsN ≈ (SN ⊗ (Hr * sparse(Diagonal(css[Nrp*Ns .+ (1:Nrp)]))))

  Asr = (QsT ⊗ Ir) * sparse(1:length(crs), 1:length(crs), crs) * (Is ⊗ Qr)
  Ars = (Is ⊗ QrT) * sparse(1:length(csr), 1:length(csr), csr) * (Qs ⊗ Ir)

  A = Arr + Ass + Ars + Asr

  Er0 = sparse([1], [1], [1], Nrp, Nrp)
  ErN = sparse([Nrp], [Nrp], [1], Nrp, Nrp)
  Es0 = sparse([1], [1], [1], Nsp, Nsp)
  EsN = sparse([Nsp], [Nsp], [1], Nsp, Nsp)

  crs0 = sparse(Diagonal(crs[1:Nrp]))
  crsN = sparse(Diagonal(crs[Nrp*Ns .+ (1:Nrp)]))
  csr0 = sparse(Diagonal(csr[1   .+ Nrp*(0:Ns)]))
  csrN = sparse(Diagonal(csr[Nrp .+ Nrp*(0:Ns)]))

  er0T = sparse([1], [1  ], [1], 1, Nrp)
  erNT = sparse([1], [Nrp], [1], 1, Nrp)
  es0T = sparse([1], [1  ], [1], 1, Nsp)
  esNT = sparse([1], [Nsp], [1], 1, Nsp)

  L1 = (Is ⊗ er0T)
  L2 = (Is ⊗ erNT)
  L3 = (es0T ⊗ Ir)
  L4 = (esNT ⊗ Ir)

  nx1 = -L1 * ys
  ny1 =  L1 * xs
  sJ1 = hypot.(nx1, ny1)
  SJ1 = Diagonal(sJ1)
  nx1 = nx1 ./ sJ1
  ny1 = ny1 ./ sJ1
  H1 = Hs
  H1I = HsI

  nx2 =  L2 * ys
  ny2 = -L2 * xs
  sJ2 = hypot.(nx2, ny2)
  SJ2 = Diagonal(sJ2)
  nx2 = nx2 ./ sJ2
  ny2 = ny2 ./ sJ2
  H2 = Hs
  H2I = HsI

  nx3 =  L3 * yr
  ny3 = -L3 * xr
  sJ3 = hypot.(nx3, ny3)
  SJ3 = Diagonal(sJ3)
  nx3 = nx3 ./ sJ3
  ny3 = ny3 ./ sJ3
  H3 = Hr
  H3I = HrI

  nx4 = -L4 * yr
  ny4 =  L4 * xr
  sJ4 = hypot.(nx4, ny4)
  SJ4 = Diagonal(sJ4)
  nx4 = nx4 ./ sJ4
  ny4 = ny4 ./ sJ4
  H4 = Hr
  H4I = HrI

  G1 =  (Is ⊗ er0T) * Sr0 + ((csr0 * Qs) ⊗ er0T)
  G2 = -(Is ⊗ erNT) * SrN - ((csrN * Qs) ⊗ erNT)
  G3 =  (es0T ⊗ Ir) * Ss0 + (es0T ⊗ (crs0 * Qr))
  G4 = -(esNT ⊗ Ir) * SsN - (esNT ⊗ (crsN * Qr))

  τsJH1 = H1 * 100 * Nsp
  τsJH2 = H2 * 100 * Nsp
  τsJH3 = H3 * 100 * Nrp
  τsJH4 = H4 * 100 * Nrp

  B1 =  (Sr0 + Sr0T) + ((csr0 * Qs + QsT * csr0) ⊗ Er0) + ((τsJH1) ⊗ Er0)
  B2 = -(SrN + SrNT) - ((csrN * Qs + QsT * csrN) ⊗ ErN) + ((τsJH2) ⊗ ErN)
  B3 =  (Ss0 + Ss0T) + (Es0 ⊗ (crs0 * Qr + QrT * crs0)) + (Es0 ⊗ (τsJH3))
  B4 = -(SsN + SsNT) - (EsN ⊗ (crsN * Qr + QrT * crsN)) + (EsN ⊗ (τsJH4))

  MIP = A + B1 + B2 + B3 + B4

  JH = sparse(1:Np, 1:Np, J) * (Hs ⊗ Hr)
  (MIP, A, (x,y), (G1, G2, G3, G4), (L1, L2, L3, L4), (H1, H2, H3, H4))
end
#}}}

#{{{ Transfinite Blend
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

function transfinite_blend(v1::T1, v2::T2, v3::T3, v4::T4, r, s;
                           e1 = (α) -> v1 * (1 .- α) / 2 + v3 * (1 .+ α) / 2,
                           e2 = (α) -> v2 * (1 .- α) / 2 + v4 * (1 .+ α) / 2,
                           e3 = (α) -> v1 * (1 .- α) / 2 + v2 * (1 .+ α) / 2,
                           e4 = (α) -> v3 * (1 .- α) / 2 + v4 * (1 .+ α) / 2
                          ) where {T1 <: Number, T2 <: Number, T3 <: Number, T4 <: Number}
  transfinite_blend(e1, e2, e3, e4, r, s)
end
#}}}

let
  p = 4
  Nr = 13
  Ns = 14

  xt = (r,s)->transfinite_blend(-1.0,  2.0, -3.0, 1.0, r, s)
  yt = (r,s)->transfinite_blend(-5.0, -1.0,  1.0, 3.0, r, s)

  (MIP, A, (x,y), (G1, G2, G3, G4), (L1, L2, L3, L4), (H1, H2, H3, H4),
   ) = locoperator(p, Nr, Ns, xt, yt)
  # (E, V) = eigen(Matrix(MIP))
  # println(extrema(E))

  τsJH1 = H1 * 100 * (Ns+1)
  τsJH2 = H2 * 100 * (Ns+1)
  τsJH3 = H3 * 100 * (Nr+1)
  τsJH4 = H4 * 100 * (Nr+1)
  S1 = G1 + τsJH1 * L1
  S2 = G2 + τsJH2 * L2
  S3 = G3 + τsJH3 * L3
  S4 = G4 + τsJH4 * L4
  M = A +
      G1' * L1 + G2' * L2 + G3' * L3 + G4' * L4 +
      L1' * S1 + L2' * S2 + L3' * S3 + L4' * S4
  @assert MIP ≈ M

  display(plot(x,y,marker=10))

  (E, V) = eigen(Matrix(M))
  println(extrema(E))

end
