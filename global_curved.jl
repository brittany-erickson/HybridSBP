if VERSION <= v"0.6.999999"
  eigen = eig
  macro isdefined(s::Symbol)
    return isdefined(s)
  end
  (!isdefined(:do_plotting)) && (do_plotting = true)
  if do_plotting
    macro plotting(ex)
      return :($(esc(ex)))
    end
  else
    macro plotting(ex)
    end
  end
  mul! = A_mul_B!
  cholesky = cholfact
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

# flatten tuples to arrays
if !@isdefined flatten_tuples
  const flatten_tuples = (x) -> reshape(collect(Iterators.flatten(x)),
                                        length(x[1]), length(x))
end

if !@isdefined ⊗
  const ⊗ = (A,B) -> kron(A, B)
end

const BC_DIRICHLET        = 1
const BC_NEUMANN          = 2
const BC_LOCKED_INTERFACE = 0
const BC_JUMP_INTERFACE   = -1

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

  τ1 = sparse(1:Nsp, 1:Nsp, 10*Nsp./sJ1)
  τ2 = sparse(1:Nsp, 1:Nsp, 10*Nsp./sJ2)
  τ3 = sparse(1:Nrp, 1:Nrp, 10*Nrp./sJ3)
  τ4 = sparse(1:Nrp, 1:Nrp, 10*Nrp./sJ4)

  # TODO: Check signs on Q terms (and update write up with correct signs)
  B1 =  (Sr0 + Sr0T) + ((csr0 * Qs + QsT * csr0) ⊗ Er0) + ((τ1 * H1 * SJ1) ⊗ Er0)
  B2 = -(SrN + SrNT) - ((csrN * Qs + QsT * csrN) ⊗ ErN) + ((τ2 * H2 * SJ2) ⊗ ErN)
  B3 =  (Ss0 + Ss0T) + (Es0 ⊗ (crs0 * Qr + QrT * crs0)) + (Es0 ⊗ (τ3 * H3 * SJ3))
  B4 = -(SsN + SsNT) - (EsN ⊗ (crsN * Qr + QrT * crsN)) + (EsN ⊗ (τ4 * H4 * SJ4))

  F1 =  (Is ⊗ er0T) * Sr0 + ((csr0 * Qs) ⊗ er0T) + ((τ1 * H1 * SJ1) ⊗ er0T)
  F2 = -(Is ⊗ erNT) * SrN - ((csrN * Qs) ⊗ erNT) + ((τ2 * H2 * SJ2) ⊗ erNT)
  F3 =  (es0T ⊗ Ir) * Ss0 + (es0T ⊗ (crs0 * Qr)) + (es0T ⊗ (τ3 * H3 * SJ3))
  F4 = -(esNT ⊗ Ir) * SsN - (esNT ⊗ (crsN * Qr)) + (esNT ⊗ (τ4 * H4 * SJ4))

  # @assert B1 ≈ F1' * L1 + L1' * F1 - ((τ1 * H1 * SJ1) ⊗ Er0)
  # @assert B2 ≈ F2' * L2 + L2' * F2 - ((τ2 * H2 * SJ2) ⊗ ErN)
  # @assert B3 ≈ F3' * L3 + L3' * F3 - (Es0 ⊗ (τ3 * H3 * SJ3))
  # @assert B4 ≈ F4' * L4 + L4' * F4 - (EsN ⊗ (τ4 * H4 * SJ4))

  M = A + B1 + B2 + B3 + B4

  # Modify the operator to handle the boundary conditions
  if !isempty(LFToB)
    F = (F1, F2, F3, F4)
    τ = (τ1, τ2, τ3, τ4)
    sJ = (sJ1, sJ2, sJ3, sJ4)
    HfI = (H1I, H2I, H3I, H4I)
    # Modify operators for the BC
    for lf = 1:4
      if LFToB[lf] == BC_NEUMANN
        M -= F[lf]' * (Diagonal(1 ./ (sJ[lf] .* diag(τ[lf]))) * HfI[lf]) * F[lf]
      elseif !(LFToB[lf] == BC_DIRICHLET ||
               LFToB[lf] == BC_LOCKED_INTERFACE ||
               LFToB[lf] == BC_JUMP_INTERFACE)
        error("invalid bc")
      end
    end
  end

  # (E, V) = eigen(Matrix(M))
  # println((minimum(E), maximum(E)))
  JH = sparse(1:Np, 1:Np, J) * (Hs ⊗ Hr)
  (M, (F1, F2, F3, F4), (L1, L2, L3, L4), (x, y), JH,
   (sJ1, sJ2, sJ3, sJ4), (nx1, nx2, nx3, nx4), (ny1, ny2, ny3, ny4),
   (H1, H2, H3, H4), (H1I, H2I, H3I, H4I), (τ1, τ2, τ3, τ4))
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
      elseif FToB[f] == BC_JUMP_INTERFACE
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
  H = sparse(1:VNp, 1:VNp, VH, VNp, VNp)

  (vstarts, M, H, X, Y, E)
end
#}}}

#{{{ gloλoperator: Build the trace operators
function gloλoperator(lop, vstarts, FToB, FToE, FToLF, EToO, EToS, Nr, Ns)
  nelems = length(lop)
  nfaces = length(FToB)
  Nλp = zeros(Int64, nfaces)
  λstarts = Array{Int64, 1}(undef, nfaces + 1)
  λstarts[1] = 1
  IT = Array{Int64,1}(undef,0)
  JT = Array{Int64,1}(undef,0)
  VT = Array{Float64,1}(undef,0)
  VD = Array{Float64,1}(undef,0)
  for f = 1:nfaces
    if FToB[f] == BC_DIRICHLET || FToB[f] == BC_NEUMANN
      λstarts[f+1] = λstarts[f]
      continue
    end
    (em, ep) = FToE[:, f]
    (fm, fp) = FToLF[:, f]
    Nλp[f] = (fm <= 2 ? Ns[em]+1 : Nr[em]+1)
    @assert Nλp[f] == (fp <= 2 ? Ns[ep]+1 : Nr[ep]+1)
    λstarts[f+1] = λstarts[f] + Nλp[f]

    @assert EToO[fm, em] && EToS[fm, em] == 1
    Fm = lop[em][2][fm]
    (Ie, Je, Ve) = findnz(Fm)
    IT = [IT; Ie .+ (λstarts[f] - 1)]
    JT = [JT; Je .+ (vstarts[em] - 1)]
    VT = [VT; Ve]

    @assert EToS[fp, ep] == 2
    Fp = lop[ep][2][fp]
    (Ie, Je, Ve) = findnz(Fp)
    # if element and face orientation do not match, then flip
    if EToO[fp, ep]
      IT = [IT; Ie .+ (λstarts[f] - 1)]
      @assert lop[em][11][fm] ≈ lop[ep][11][fp]
      @assert lop[em][6][fm] ≈ lop[ep][6][fp]
      @assert lop[em][9][fm] ≈ lop[ep][9][fp]
    else
      IT = [IT; λstarts[f+1] .- Ie]
      @assert lop[em][11][fm] ≈ rot180(lop[ep][11][fp])
      @assert lop[em][6][fm] ≈ lop[ep][6][fp][end:-1:1]
      @assert lop[em][9][fm] ≈ rot180(lop[ep][9][fp])
    end
    JT = [JT; Je .+ (vstarts[ep] - 1)]
    VT = [VT; Ve]

    sJf = lop[em][6][fm]
    Hf = Vector(diag(lop[em][9][fm]))
    τf = Vector(diag(lop[em][11][fm]))
    VD = [VD; 2 * sJf .* Hf .* τf]

  end
  λNp = λstarts[nfaces+1]-1
  VNp = vstarts[nelems+1]-1
  T = sparse(IT, JT, VT, λNp, VNp)
  # Ttranspose = sparse(JT, IT, VT, VNp, λNp)
  D = sparse(Diagonal(VD))
  (λstarts, T, D)
end
#}}}

#{{{ volbcarray()
function locbcarray!(ge, lop, LFToB, bc_Dirichlet, bc_Neumann, in_jump,
                     bcargs = ())
  (~, F, L, (x, y), ~, sJ, nx, ny, ~, ~, τ) = lop
  ge[:] .= 0
  for lf = 1:4
    (xf, yf) = (L[lf] * x, L[lf] * y)
    if LFToB[lf] == BC_DIRICHLET
      vf = bc_Dirichlet(lf, xf, yf, bcargs...)
    elseif LFToB[lf] == BC_NEUMANN
      vf = bc_Neumann(lf, xf, yf, nx[lf], ny[lf], bcargs...) ./ diag(τ[lf])
    elseif LFToB[lf] == BC_LOCKED_INTERFACE
      continue # nothing to do here
    elseif LFToB[lf] == BC_JUMP_INTERFACE
      # In this case we need to add in half the jump
      vf = in_jump(lf, xf, yf, bcargs...) / 2
    else
      error("invalid bc")
    end
    ge[:] += F[lf]' * vf
  end
end
#}}}

#{{{ cg
function enorm(A::UniformScaling{T}, g, tmp) where T
  g' * g
end
function enorm(M::Vector{T}, g, tmp) where T
    @. tmp = M * g
    g' * tmp
end
function cg(u0, b, A; tol=1e-8, MaxIter=100, M = I)
  cg(u0, b, (y,x)->mul!(y, A, x); tol=tol, MaxIter=MaxIter, M=M)
end
function cg(u0, b, A::Function; tol=1e-8, MaxIter=100, M = I)
  u = copy(u0)
  w = similar(u0)
  d = similar(u0)
  g = similar(u0)
  tmp = similar(u0)
  k = cg!(u, w, d, g, tmp, b, A; tol=tol, MaxIter=MaxIter, M=M)
  (u, k)
end
function cg!(u, w, d, g, tmp, b, A::Function; tol=1e-8, MaxIter=100, M = I)

  A(w, u)
  @. d = b - w
  @. g = -d

  gkTgk = g' * g

  err = enorm(M, g, tmp)
  nmx = enorm(M, u, tmp)
  # err = g' * g
  # nmx = u' * u
  tol2 = tol^2
  for k = 1:MaxIter
    A(w, d)

    alpha = gkTgk / (d' * w)

    @. u = u + alpha * d

    @. g = g + alpha * w

    gk1Tgk1 = g' * g

    beta = gk1Tgk1 / gkTgk

    @. d = -g + beta * d

    gkTgk = gk1Tgk1

    err = enorm(M, g, tmp)
    nmx = enorm(M, u, tmp)
    # err = g' * g
    # nmx = u' * u
    if err < tol2 * (1 + nmx)
      return k
    end
  end
  -MaxIter
end
function testcg(N)
  (Q,R) = qr(rand(N,N))
  A = Q * Diagonal(rand(N)) * Q'
  A = (A+A')/2
  x = rand(N)

  b = A*x
  f = x->A*x
  b = A*x

  u = A*x
  (x0,k) = cg(u,b,A, tol=1e-6, MaxIter=100, M=I)
  println((norm(x0 - x), k))

  u = A*x
  (x0,k) = cg(u,b,A, tol=1e-6, MaxIter=100, M=I)
  println((norm(x0 - x), k))

  nothing
end
#}}}
