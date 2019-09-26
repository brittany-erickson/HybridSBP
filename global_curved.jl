using SparseArrays
using LinearAlgebra

if do_plotting
  macro plotting(ex)
    return :($(esc(ex)))
  end
else
  macro plotting(ex)
  end
end

@plotting let
  using Plots
  # pyplot()
end

include("diagonal_sbp.jl")

# flatten tuples to arrays
flatten_tuples(x) = reshape(collect(Iterators.flatten(x)), length(x[1]),
                            length(x))

⊗(A,B) = kron(A, B)

const BC_DIRICHLET        = 1
const BC_NEUMANN          = 2
const BC_LOCKED_INTERFACE = 0
const BC_JUMP_INTERFACE   = 7

#{{{ Transfinite Blend
function transfinite_blend(α1, α2, α3, α4, r, s)
  # +---4---+
  # |       |
  # 1       2
  # |       |
  # +---3---+
  @assert [α1(-1) α2(-1) α1( 1) α2( 1)] ≈ [α3(-1) α3( 1) α4(-1) α4( 1)]

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

#{{{ connectivityarrays
function connectivityarrays(EToV, EToF)
  # number of elements
  nelems = size(EToV, 2)
  nfaces = maximum(maximum(EToF))

  # Determine secondary arrays
  # FToE : Unique Global Face to Element Number
  # FToLF: Unique Global Face to Element local face number
  # EToO : Element to Unique Global Faces Orientation
  # EToS : Element to Unique Global Face Side

  FToE  = zeros(Int64, 2, nfaces)
  FToLF = zeros(Int64, 2, nfaces)
  EToO  = Array{Bool,2}(undef, 4, nelems)
  EToS  = zeros(Int64, 4, nelems)

  # Local Face to Local Vertex map
  LFToLV = flatten_tuples(((1,3), (2, 4), (1,2), (3,4)))
  for e = 1:nelems
    for lf = 1:4
      gf = EToF[lf, e]
      if FToE[1, gf] == 0
        @assert FToLF[1, gf] == 0
        FToE[1, gf] = e
        FToLF[1, gf] = lf
        EToO[lf, e] = true
        EToS[lf, e] = 1
      else
        @assert FToE[2, gf] == 0
        @assert FToLF[2, gf] == 0
        FToE[2, gf] = e
        FToLF[2, gf] = lf
        EToS[lf, e] = 2

        ne = FToE[1, gf]
        nf = FToLF[1, gf]

        nv = EToV[LFToLV[:,nf], ne]
        lv = EToV[LFToLV[:,lf], e]
        if nv == lv
          EToO[lf, e] = true
        elseif nv[end:-1:1] == lv
          EToO[lf, e] = false
        else
          error("problem with connectivity")
        end
      end
    end
  end
  (FToE, FToLF, EToO, EToS)
end
#}}}

#{{{ locoperator
function locoperator(p, Nr, Ns, xf, yf; pm = p+2, LFToB = [], τscale = 2)
  Nrp = Nr + 1
  Nsp = Ns + 1
  Np = Nrp * Nsp

  # Derivative operators for the metric terms
  (DrM, _, _, _) = diagonal_sbp_D1(pm, Nr; xc = (-1,1))
  (DsM, _, _, _) = diagonal_sbp_D1(pm, Ns; xc = (-1,1))

  # Derivative operators for the rest of the computation
  (Dr, HrI, Hr, r) = diagonal_sbp_D1(p, Nr; xc = (-1,1))
  Qr = Hr * Dr
  QrT = sparse(transpose(Qr))

  (Ds, HsI, Hs, s) = diagonal_sbp_D1(p, Ns; xc = (-1,1))
  Qs = Hs * Ds
  QsT = sparse(transpose(Qs))

  # Identity matrices for the comuptation
  Ir = sparse(I, Nrp, Nrp)
  Is = sparse(I, Nsp, Nsp)

  # Create the mesh
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

  # variable coefficient matrix components
  crr = J .* (rx .* rx + ry .* ry)
  crs = csr = J .* (sx .* rx + sy .* ry)
  css = J .* (sx .* sx + sy .* sy)

  #{{{ Set up the rr derivative matrix
  ISr0 = Array{Int64,1}(undef,0)
  JSr0 = Array{Int64,1}(undef,0)
  VSr0 = Array{Float64,1}(undef,0)
  ISrN = Array{Int64,1}(undef,0)
  JSrN = Array{Int64,1}(undef,0)
  VSrN = Array{Float64,1}(undef,0)

  (_, S0e, SNe, _, _, Ae, _) = variable_diagonal_sbp_D2(p, Nr, rand(Nrp))
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
    (_, S0e, SNe, _, _, Ae, _) = variable_diagonal_sbp_D2(p, Nr, crr[rng])
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
  Ãrr = sparse(IArr[1:stArr], JArr[1:stArr], VArr[1:stArr], Np, Np)
  Sr0 = sparse(ISr0[1:stSr0], JSr0[1:stSr0], VSr0[1:stSr0], Np, Np)
  SrN = sparse(ISrN[1:stSrN], JSrN[1:stSrN], VSrN[1:stSrN], Np, Np)
  Sr0T = sparse(JSr0[1:stSr0], ISr0[1:stSr0], VSr0[1:stSr0], Np, Np)
  SrNT = sparse(JSrN[1:stSrN], ISrN[1:stSrN], VSrN[1:stSrN], Np, Np)
  #= affine mesh test
  # @assert Ãrr ≈ Ãrr'
  (D2, S0, SN, _, _, _) = diagonal_sbp_D2(p, Nr)
  Ar = SN - S0 - Hr * D2
  @assert Ãrr ≈ Hs ⊗ Ar
  =#
  # @assert Sr0 ≈ ((sparse(Diagonal(crr[1   .+ Nrp*(0:Ns)])) * Hs) ⊗ S0)
  # @assert SrN ≈ ((sparse(Diagonal(crr[Nrp .+ Nrp*(0:Ns)])) * Hs) ⊗ SN)
  #}}}

  #{{{ Set up the ss derivative matrix
  (_, S0e, SNe, _, _, Ae, _) = variable_diagonal_sbp_D2(p, Ns, rand(Nsp))
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
    (_, S0e, SNe, _, _, Ae, _) = variable_diagonal_sbp_D2(p, Ns, css[rng])

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
  Ãss = sparse(IAss[1:stAss], JAss[1:stAss], VAss[1:stAss], Np, Np)
  Ss0 = sparse(ISs0[1:stSs0], JSs0[1:stSs0], VSs0[1:stSs0], Np, Np)
  SsN = sparse(ISsN[1:stSsN], JSsN[1:stSsN], VSsN[1:stSsN], Np, Np)
  Ss0T = sparse(JSs0[1:stSs0], ISs0[1:stSs0], VSs0[1:stSs0], Np, Np)
  SsNT = sparse(JSsN[1:stSsN], ISsN[1:stSsN], VSsN[1:stSsN], Np, Np)
  # @assert Ãss ≈ Ãss'
  #= affine mesh test
  (D2, S0, SN, _, _, _) = diagonal_sbp_D2(p, Ns)
  As = SN - S0 - Hs * D2
  @assert Ãss ≈ As ⊗ Hr
  =#
  # @assert Ss0 ≈ (S0 ⊗ (Hr * sparse(Diagonal(css[1:Nrp]))))
  # @assert SsN ≈ (SN ⊗ (Hr * sparse(Diagonal(css[Nrp*Ns .+ (1:Nrp)]))))
  #}}}

  #{{{ Set up the sr and rs derivative matrices
  Ãsr = (QsT ⊗ Ir) * sparse(1:length(crs), 1:length(crs), crs) * (Is ⊗ Qr)
  Ãrs = (Is ⊗ QrT) * sparse(1:length(csr), 1:length(csr), csr) * (Qs ⊗ Ir)
  #}}}

  Ã = Ãrr + Ãss + Ãrs + Ãsr

  #
  # Boundary point matrices
  #
  Er0 = sparse([1], [1], [1], Nrp, Nrp)
  ErN = sparse([Nrp], [Nrp], [1], Nrp, Nrp)
  Es0 = sparse([1], [1], [1], Nsp, Nsp)
  EsN = sparse([Nsp], [Nsp], [1], Nsp, Nsp)

  er0T = sparse([1], [1  ], [1], 1, Nrp)
  erNT = sparse([1], [Nrp], [1], 1, Nrp)
  es0T = sparse([1], [1  ], [1], 1, Nsp)
  esNT = sparse([1], [Nsp], [1], 1, Nsp)

  L1 = (Is ⊗ er0T)
  L2 = (Is ⊗ erNT)
  L3 = (es0T ⊗ Ir)
  L4 = (esNT ⊗ Ir)

  #
  # Store coefficient matrices as matrices
  #
  crs0 = sparse(Diagonal(crs[1:Nrp]))
  crsN = sparse(Diagonal(crs[Nrp*Ns .+ (1:Nrp)]))
  csr0 = sparse(Diagonal(csr[1   .+ Nrp*(0:Ns)]))
  csrN = sparse(Diagonal(csr[Nrp .+ Nrp*(0:Ns)]))

  #
  # Block surface matrices
  #
  nx1 = -L1 * ys
  ny1 =  L1 * xs
  sJ1 = hypot.(nx1, ny1)
  nx1 = nx1 ./ sJ1
  ny1 = ny1 ./ sJ1
  H1 = Hs
  H1I = HsI

  nx2 =  L2 * ys
  ny2 = -L2 * xs
  sJ2 = hypot.(nx2, ny2)
  nx2 = nx2 ./ sJ2
  ny2 = ny2 ./ sJ2
  H2 = Hs
  H2I = HsI

  nx3 =  L3 * yr
  ny3 = -L3 * xr
  sJ3 = hypot.(nx3, ny3)
  nx3 = nx3 ./ sJ3
  ny3 = ny3 ./ sJ3
  H3 = Hr
  H3I = HrI

  nx4 = -L4 * yr
  ny4 =  L4 * xr
  sJ4 = hypot.(nx4, ny4)
  nx4 = nx4 ./ sJ4
  ny4 = ny4 ./ sJ4
  H4 = Hr
  H4I = HrI

  #
  # Penalty terms
  #
  if p == 2
    l = 2
    β = 0.363636363
    α = 1 / 2
  elseif p == 4
    l = 4
    β = 0.2505765857
    α = 17 / 48
  elseif p == 6
    l = 6
    β = 0.1878687080
    α = 13649 / 43200
  else
    error("unknown order")
  end

  crr = reshape(crr, Nrp, Nsp)
  css = reshape(css, Nrp, Nsp)
  crs = reshape(crs, Nrp, Nsp)
  csr = reshape(csr, Nrp, Nsp)
  ψmin = reshape((crr + css - sqrt.((crr - css).^2 + 4crs.^2)) / 2, Nrp, Nsp)

  hr = 2 / Nr
  hs = 2 / Ns

  ψ1 = ψmin[  1, :]
  ψ2 = ψmin[Nrp, :]
  ψ3 = ψmin[:,   1]
  ψ4 = ψmin[:, Nsp]
  for k = 2:l
    ψ1 = min.(ψ1, ψmin[k, :])
    ψ2 = min.(ψ2, ψmin[Nr-k, :])
    ψ3 = min.(ψ3, ψmin[:, k])
    ψ4 = min.(ψ4, ψmin[:, Ns-k])
  end
  τ1 = (2τscale / hr) * (crr[  1, :].^2 / β + crs[  1, :].^2 / α) ./ ψ1
  τ2 = (2τscale / hr) * (crr[Nrp, :].^2 / β + crs[Nrp, :].^2 / α) ./ ψ2
  τ3 = (2τscale / hs) * (css[:,   1].^2 / β + crs[:,   1].^2 / α) ./ ψ3
  τ4 = (2τscale / hs) * (css[:, Nsp].^2 / β + crs[:, Nsp].^2 / α) ./ ψ4

  τ1 = sparse(1:Nsp, 1:Nsp, τ1)
  τ2 = sparse(1:Nsp, 1:Nsp, τ2)
  τ3 = sparse(1:Nrp, 1:Nrp, τ3)
  τ4 = sparse(1:Nrp, 1:Nrp, τ4)

  C̃1 =  (Sr0 + Sr0T) + ((csr0 * Qs + QsT * csr0) ⊗ Er0) + ((τ1 * H1) ⊗ Er0)
  C̃2 = -(SrN + SrNT) - ((csrN * Qs + QsT * csrN) ⊗ ErN) + ((τ2 * H2) ⊗ ErN)
  C̃3 =  (Ss0 + Ss0T) + (Es0 ⊗ (crs0 * Qr + QrT * crs0)) + (Es0 ⊗ (τ3 * H3))
  C̃4 = -(SsN + SsNT) - (EsN ⊗ (crsN * Qr + QrT * crsN)) + (EsN ⊗ (τ4 * H4))

  F1 =  (Is ⊗ er0T) * Sr0 + ((csr0 * Qs) ⊗ er0T) + ((τ1 * H1) ⊗ er0T)
  F2 = -(Is ⊗ erNT) * SrN - ((csrN * Qs) ⊗ erNT) + ((τ2 * H2) ⊗ erNT)
  F3 =  (es0T ⊗ Ir) * Ss0 + (es0T ⊗ (crs0 * Qr)) + (es0T ⊗ (τ3 * H3))
  F4 = -(esNT ⊗ Ir) * SsN - (esNT ⊗ (crsN * Qr)) + (esNT ⊗ (τ4 * H4))

  G1 =  (Is ⊗ er0T) * Sr0 + ((csr0 * Qs) ⊗ er0T)
  G2 = -(Is ⊗ erNT) * SrN - ((csrN * Qs) ⊗ erNT)
  G3 =  (es0T ⊗ Ir) * Ss0 + (es0T ⊗ (crs0 * Qr))
  G4 = -(esNT ⊗ Ir) * SsN - (esNT ⊗ (crsN * Qr))

  # @assert C̃1 ≈ F1' * L1 + L1' * F1 - ((τ1 * H1) ⊗ Er0)
  # @assert C̃2 ≈ F2' * L2 + L2' * F2 - ((τ2 * H2) ⊗ ErN)
  # @assert C̃3 ≈ F3' * L3 + L3' * F3 - (Es0 ⊗ (τ3 * H3))
  # @assert C̃4 ≈ F4' * L4 + L4' * F4 - (EsN ⊗ (τ4 * H4))

  M̃ = Ã + C̃1 + C̃2 + C̃3 + C̃4

  # Modify the operator to handle the boundary conditions
  if !isempty(LFToB)
    F = (F1, F2, F3, F4)
    τ = (τ1, τ2, τ3, τ4)
    HfI = (H1I, H2I, H3I, H4I)
    # Modify operators for the BC
    for lf = 1:4
      if LFToB[lf] == BC_NEUMANN
        M̃ -= F[lf]' * (Diagonal(1 ./ (diag(τ[lf]))) * HfI[lf]) * F[lf]
      elseif !(LFToB[lf] == BC_DIRICHLET ||
               LFToB[lf] == BC_LOCKED_INTERFACE ||
               LFToB[lf] >= BC_JUMP_INTERFACE)
        error("invalid bc")
      end
    end
  end

  # (E, V) = eigen(Matrix(M̃))
  # println((minimum(E), maximum(E)))
  JH = sparse(1:Np, 1:Np, J) * (Hs ⊗ Hr)
  (M̃, (F1, F2, F3, F4), (L1, L2, L3, L4), (x, y), JH,
   (sJ1, sJ2, sJ3, sJ4), (nx1, nx2, nx3, nx4), (ny1, ny2, ny3, ny4),
   (H1, H2, H3, H4), (H1I, H2I, H3I, H4I), (τ1, τ2, τ3, τ4), (G1, G2, G3, G4))
end
#}}}

#{{{ plot the mesh
@plotting function plotmesh(p, lop, Nr, Ns, EToF, FToB)
  for e = 1:length(lop)
    (x, y) = lop[e][4]
    plot!(p, reshape(x, Nr[e]+1, Ns[e]+1), reshape(y, Nr[e]+1, Ns[e]+1),
          color=:black, legend=:none)
    plot!(p, reshape(x, Nr[e]+1, Ns[e]+1)', reshape(y, Nr[e]+1, Ns[e]+1)',
          color=:black, legend=:none)
    L = lop[e][3]
    for lf = 1:4
      f = EToF[lf, e]
      if FToB[f] == BC_DIRICHLET
        plot!(p, L[lf] * x, L[lf] * y, color=:red, legend=:none, linewidth=3)
      elseif FToB[f] == BC_NEUMANN
        plot!(p, L[lf] * x, L[lf] * y, color=:blue, legend=:none, linewidth=3)
      elseif FToB[f] == BC_LOCKED_INTERFACE
        plot!(p, L[lf] * x, L[lf] * y, color=:green, legend=:none, linewidth=3)
      elseif FToB[f] >= BC_JUMP_INTERFACE
        plot!(p, L[lf] * x, L[lf] * y, color=:purple, legend=:none, linewidth=3)
      else
        error("invalid bc")
      end
    end
  end
end
#}}}

#{{{glovoloperator: Assemble the global volume operators
function glovoloperator(lop, Nr, Ns)
  nelems = length(lop)
  vstarts = Array{Int64, 1}(undef, nelems + 1)
  vstarts[1] = 1
  Np = Array{Int64, 1}(undef, nelems)
  IM = Array{Int64,1}(undef,0)
  JM = Array{Int64,1}(undef,0)
  VM = Array{Float64,1}(undef,0)
  VH = Array{Float64,1}(undef,0)
  X = Array{Float64,1}(undef,0)
  Y = Array{Float64,1}(undef,0)
  E = Array{Float64,1}(undef,0)
  for e = 1:nelems
    # Fill arrays to build global sparse matrix
    Np[e] = (Nr[e]+1)*(Ns[e]+1)
    vstarts[e+1] = vstarts[e] + Np[e]
    M = lop[e][1]
    (Ie, Je, Ve) = findnz(M)
    IM = [IM;Ie .+ (vstarts[e]-1)]
    JM = [JM;Je .+ (vstarts[e]-1)]
    VM = [VM;Ve]

    # Global "mass" matrix
    H = lop[e][5]
    VH = [VH;Vector(diag(H))]

    # global coordinates and element number array (needed for jump)
    x = lop[e][4][1]
    X = [X;x]
    y = lop[e][4][2]
    Y = [Y;y]
    E = [E;e * ones(Np[e])]

  end
  VNp = vstarts[nelems+1]-1 # total number of volume points
  M = sparse(IM, JM, VM, VNp, VNp)

  (vstarts, M, VH, X, Y, E)
end
#}}}

#{{{ gloλoperator: Build the trace operators
function gloλoperator(lop, vstarts, FToB, FToE, FToLF, EToO, EToS, Nr, Ns)
  nelems = length(lop)
  nfaces = length(FToB)
  Nλp = zeros(Int64, nfaces)
  FToλstarts = Array{Int64, 1}(undef, nfaces + 1)
  FToλstarts[1] = 1
  IT = Array{Int64,1}(undef,0)
  JT = Array{Int64,1}(undef,0)
  VT = Array{Float64,1}(undef,0)
  VD = Array{Float64,1}(undef,0)
  for f = 1:nfaces
    if FToB[f] == BC_DIRICHLET || FToB[f] == BC_NEUMANN
      FToλstarts[f+1] = FToλstarts[f]
      continue
    end
    (em, ep) = FToE[:, f]
    (fm, fp) = FToLF[:, f]
    Nλp[f] = (fm <= 2 ? Ns[em]+1 : Nr[em]+1)
    @assert Nλp[f] == (fp <= 2 ? Ns[ep]+1 : Nr[ep]+1)
    FToλstarts[f+1] = FToλstarts[f] + Nλp[f]

    @assert EToO[fm, em] && EToS[fm, em] == 1
    Fm = lop[em][2][fm]
    (Ie, Je, Ve) = findnz(Fm)
    IT = [IT; Ie .+ (FToλstarts[f] - 1)]
    JT = [JT; Je .+ (vstarts[em] - 1)]
    VT = [VT; Ve]

    @assert EToS[fp, ep] == 2
    Fp = lop[ep][2][fp]
    (Ie, Je, Ve) = findnz(Fp)
    # if element and face orientation do not match, then flip
    if EToO[fp, ep]
      IT = [IT; Ie .+ (FToλstarts[f] - 1)]
      # @assert lop[em][11][fm] ≈ lop[ep][11][fp]
      # @assert lop[em][6][fm] ≈ lop[ep][6][fp]
      # @assert lop[em][9][fm] ≈ lop[ep][9][fp]
      τm = Vector(diag(lop[em][11][fm]))
      τp = Vector(diag(lop[ep][11][fp]))
    else
      IT = [IT; FToλstarts[f+1] .- Ie]
      # @assert lop[em][11][fm] ≈ rot180(lop[ep][11][fp])
      # @assert lop[em][6][fm] ≈ lop[ep][6][fp][end:-1:1]
      # @assert lop[em][9][fm] ≈ rot180(lop[ep][9][fp])
      τm = Vector(diag(lop[em][11][fm]))
      τp = Vector(diag(rot180(lop[ep][11][fp])))
    end
    JT = [JT; Je .+ (vstarts[ep] - 1)]
    VT = [VT; Ve]

    Hf = Vector(diag(lop[em][9][fm]))
    VD = [VD; Hf .* (τm + τp)]

  end
  λNp = FToλstarts[nfaces+1]-1
  VNp = vstarts[nelems+1]-1
  T = sparse(IT, JT, VT, λNp, VNp)
  # Ttranspose = sparse(JT, IT, VT, VNp, λNp)
  (FToλstarts, T, VD)
end
#}}}

#{{{ volbcarray()
function locbcarray!(ge, lop, LFToB, bc_Dirichlet, bc_Neumann, in_jump,
                     bcargs = ())
  (_, F, L, (x, y), _, sJ, nx, ny, _, _, τ) = lop
  ge[:] .= 0
  for lf = 1:4
    (xf, yf) = (L[lf] * x, L[lf] * y)
    if LFToB[lf] == BC_DIRICHLET
      vf = bc_Dirichlet(lf, xf, yf, bcargs...)
    elseif LFToB[lf] == BC_NEUMANN
      gN = bc_Neumann(lf, xf, yf, nx[lf], ny[lf], bcargs...)
      vf = sJ[lf] .* gN ./ diag(τ[lf])
    elseif LFToB[lf] == BC_LOCKED_INTERFACE
      continue # nothing to do here
    elseif LFToB[lf] >= BC_JUMP_INTERFACE
      # In this case we need to add in half the jump
      vf = in_jump(lf, xf, yf, bcargs...) / 2
    else
      error("invalid bc")
    end
    ge[:] += F[lf]' * vf
  end
end
#}}}

#{{{
struct SBPLocalOperator1{T<:Real, S<:Factorization}
  offset::Array{Int64,1}
  H::Array{T,1}
  X::Array{T,1}
  Y::Array{T,1}
  E::Array{Int64,1}
  F::Array{S,1}
  SBPLocalOperator1{T,S}(vstarts::Array{Int64,1}, H::Array{T,1}, X::Array{T,1},
                         Y::Array{T,1}, E::Array{Int64,1},
                         F::Array{S,1}) where {T<:Real, S<:Factorization} =
  new(vstarts, H, X, Y, E, F)
end

function SBPLocalOperator1(lop, Nr, Ns, factorization)
  nelems = length(lop)
  vstarts = Array{Int64, 1}(undef, nelems + 1)
  vstarts[1] = 1
  Np = Array{Int64, 1}(undef, nelems)
  VH = Array{Float64,1}(undef,0)
  X = Array{Float64,1}(undef,0)
  Y = Array{Float64,1}(undef,0)
  E = Array{Int64,1}(undef,0)
  FTYPE = typeof(factorization(sparse([1],[1],[1.0])))
  factors = Array{FTYPE, 1}(undef, nelems)
  for e = 1:nelems
    # Fill arrays to build global sparse matrix
    Np[e] = (Nr[e]+1)*(Ns[e]+1)
    vstarts[e+1] = vstarts[e] + Np[e]

    # Global "mass" matrix
    H = lop[e][5]
    VH = [VH;Vector(diag(H))]

    # global coordinates and element number array (needed for jump)
    x = lop[e][4][1]
    X = [X;x]
    y = lop[e][4][2]
    Y = [Y;y]
    E = [E;e * ones(Int64, Np[e])]

    factors[e] = factorization(lop[e][1])
  end
  VNp = vstarts[nelems+1]-1 # total number of volume points

  SBPLocalOperator1{Float64, FTYPE}(vstarts, VH, X, Y, E, factors)
end
#}}}

function LocalGlobalOperators(lop, Nr, Ns, FToB, FToE, FToLF, EToO, EToS,
                              factorization)
  M = SBPLocalOperator1(lop, Nr, Ns, factorization)
  (FToλstarts, FbarT, D) = gloλoperator(lop, M.offset, FToB, FToE, FToLF, EToO,
                                        EToS, Nr, Ns)
  (M, FbarT, D, M.offset, FToλstarts)
end

function bcstarts(FToB, FToE, FToLF, bc_type, Nr, Ns)
  nfaces = length(FToB)
  bcstarts = Array{Int64, 1}(undef, nfaces + 1)
  bcstarts[1] = 1
  for f = 1:nfaces
    if FToB[f] ∈ bc_type
      e  = FToE[1,f]
      lf = FToLF[1,f]
      bcstarts[f+1] = bcstarts[f] + (lf ∈ (1,2) ? Ns[e] : Nr[e]) + 1
    else
      bcstarts[f+1] = bcstarts[f]
    end
  end
  bcstarts
end

function LocalToGLobalRHS!(b, g, u, M, FbarT, vstarts, lockedblock)
  for e = 1:length(M)
    if !lockedblock[e]
      @views u[vstarts[e]:(vstarts[e+1]-1)] = M[e] \ g[vstarts[e]:(vstarts[e+1]-1)]
      #=
      ldiv!((@view u[vstarts[e]:(vstarts[e+1]-1)]), M[e],
            (@view g[vstarts[e]:(vstarts[e+1]-1)]))
      =#
    else
      @views u[vstarts[e]:(vstarts[e+1]-1)] .= 0
    end
  end
  mul!(b, FbarT, u)
end

#{{{ assembleλmatrix: Schur complement system
function assembleλmatrix(FToλstarts, vstarts, EToF, FToB, F, D, FbarT)
  nfaces = length(FToλstarts)-1
  nelems = length(vstarts)-1
  λNp = FToλstarts[nfaces+1]-1
  sz = λNp

  for e = 1:nelems
    lλs = Array{Int64, 1}(undef, 4)
    for lf = 1:4
      f = EToF[lf,e]
      lλs[lf] = FToλstarts[f+1] - FToλstarts[f]
    end
    for lf = 1:4
      sz += lλs[lf]*sum(lλs)
    end
  end
  Ie = Array{Int64, 1}(undef, sz)
  Je = Array{Int64, 1}(undef, sz)
  Ve = Array{Float64, 1}(undef, sz)
  Ie[1:λNp] = 1:λNp
  Je[1:λNp] = 1:λNp
  Ve[1:λNp] = D
  offset = λNp
  Fbar = FbarT'
  for e = 1:nelems
    # println((e, nelems))
    vrng = vstarts[e]:(vstarts[e+1]-1)
    for lf = 1:4
      f = EToF[lf,e]
      if FToB[f] == BC_LOCKED_INTERFACE || FToB[f] >= BC_JUMP_INTERFACE
        λrng = FToλstarts[f]:(FToλstarts[f+1]-1)
        B = Matrix(F[e] \ Fbar[vrng, λrng])
        for lf2 = 1:4
          f2 = EToF[lf2,e]
          if FToB[f2] == BC_LOCKED_INTERFACE || FToB[f2] >= BC_JUMP_INTERFACE
            λrng2 = FToλstarts[f2]:(FToλstarts[f2+1]-1)
            C = FbarT[λrng2, vrng] * B
            λblck = λrng*ones(Int64, 1, length(λrng2))
            λblck2 = ones(Int64, length(λrng), 1) * λrng2'
            last = length(λrng) * length(λrng2)
            Ie[offset.+(1:last)] = λblck[:]
            Je[offset.+(1:last)] = λblck2[:]
            Ve[offset.+(1:last)] = -C'[:]
            offset += last
          end
        end
      end
    end
  end
  @assert offset == sz
  B = sparse(Ie, Je, Ve, λNp, λNp)
  @assert B ≈ B'
  # println((λNp * λNp, nnz(B), nnz(B) / λNp^2))
  B
end

#}}}

# {{{ Constructor for inp files
function read_inp_2d(T, S, filename::String; bc_map=1:10000)
  # {{{ Read in the file
  f = try
    open(filename)
  catch
    error("InpRead cannot open \"$filename\" ")
  end
  lines = readlines(f)
  close(f)
  # }}}

  # {{{ Read in nodes
  str = "NSET=ALLNODES"
  linenum = SeekToSubstring(lines, str);
  linenum > 0 || error("did not find: $str")
  num_nodes = 0
  for l = linenum+1:length(lines)
    occursin(r"^\s*[0-9]*\s*,.*", lines[l]) ? num_nodes+=1 : break
  end
  Vx = fill(S(NaN), num_nodes)
  Vy = fill(S(NaN), num_nodes)
  Vz = fill(S(NaN), num_nodes)
  for l = linenum .+ (1:num_nodes)
    node_data = split(lines[l], r"\s|,", keepempty=false)
    (node_num, node_x, node_y, node_z) = try
      (parse(T, node_data[1]),
       parse(S, node_data[2]),
       parse(S, node_data[3]),
       parse(S, node_data[4]))
    catch
      error("cannot parse line $l: \"$(lines[l])\" ")
    end

    Vx[node_num] = node_x
    Vy[node_num] = node_y
    Vz[node_num] = node_z
  end
  # }}}

  # {{{ Read in Elements
  str = "ELEMENT"
  linenum = SeekToSubstring(lines, str);
  num_elm = 0
  while linenum > 0
    for l = linenum .+ (1:length(lines))
      occursin(r"^\s*[0-9]*\s*,.*", lines[l]) ? num_elm+=1 : break
    end
    linenum = SeekToSubstring(lines, str; first=linenum+1)
  end
  num_elm > 0 || error("did not find any element")

  EToV = fill(T(0), 4, num_elm)
  EToBlock = fill(T(0), num_elm)
  linenum = SeekToSubstring(lines, str);
  while linenum > 0
    foo = split(lines[linenum], r"[^0-9]", keepempty=false)
    B = parse(T, foo[end])
    for l = linenum .+ (1:num_elm)
      elm_data = split(lines[l], r"\s|,", keepempty=false)
      # read into z-order
      (elm_num, elm_v1, elm_v2, elm_v4, elm_v3) = try
        (parse(T, elm_data[1]),
         parse(T, elm_data[2]),
        parse(T, elm_data[3]),
        parse(T, elm_data[4]),
        parse(T, elm_data[5]))
      catch
        break
      end
      EToV[:, elm_num] = [elm_v1, elm_v2, elm_v3, elm_v4]
      EToBlock[elm_num] = B
    end
    linenum = SeekToSubstring(lines, str; first=linenum+1)
  end
  # }}}

  # {{{ Determine connectivity
  EToF = fill(T(0), 4, num_elm)

  VsToF = Dict{Tuple{Int64, Int64}, Int64}()
  numfaces = 0
  for e = 1:num_elm
    for lf = 1:4
      if lf == 1
        Vs = (EToV[1, e], EToV[3, e])
      elseif lf == 2
        Vs = (EToV[2, e], EToV[4, e])
      elseif lf == 3
        Vs = (EToV[1, e], EToV[2, e])
      elseif lf == 4
        Vs = (EToV[3, e], EToV[4, e])
      end
      if Vs[1] > Vs[2]
        Vs = (Vs[2], Vs[1])
      end
      if haskey(VsToF, Vs)
        EToF[lf, e] = VsToF[Vs]
      else
        numfaces = numfaces + 1
        EToF[lf, e] = VsToF[Vs] = numfaces
      end
    end
  end
  #}}}

  # {{{ Read in side set info
  FToB = Array{T, 1}(undef, numfaces)
  fill!(FToB, BC_LOCKED_INTERFACE)
  linenum = SeekToSubstring(lines, "\\*ELSET")
  inp_to_zorder = [3,  2, 4, 1]
  while linenum > 0
    foo = split(lines[linenum], r"[^0-9]", keepempty=false)
    (bc, face) = try
      (parse(T, foo[1]),
       parse(T, foo[2]))
    catch
      error("cannot parse line $linenum: \"$(lines[linenum])\" ")
    end
    bc = bc_map[bc]
    face = inp_to_zorder[face]
    for l = linenum+1:length(lines)
      if !occursin(r"^\s*[0-9]+", lines[l])
        break
      end
      elms = split(lines[l], r"\s|,", keepempty=false)
      for elm in elms
        elm = try
          parse(T, elm)
        catch
          error("cannot parse line $linenum: \"$(lines[l])\" ")
        end
        if bc == 3
          bc = BC_LOCKED_INTERFACE
        end
        FToB[EToF[face, elm]] = bc
        @assert (bc == BC_DIRICHLET || bc == BC_NEUMANN ||
                 bc == BC_LOCKED_INTERFACE || bc >= BC_JUMP_INTERFACE)
      end
    end
    linenum = SeekToSubstring(lines, "\\*ELSET"; first=linenum+1)
  end
  # }}}

  ([Vx Vy]', EToV, EToF, FToB, EToBlock)
end
read_inp_2d(filename;kw...) = read_inp_2d(Int64, Float64, filename;kw...)

function SeekToSubstring(lines, substring; first=1)
  for l = first:length(lines)
    if occursin(Regex(".*$(substring).*"), lines[l])
      return l
    end
  end
  return -1
end

# }}}
