if VERSION == v"0.6.3"
  eigen = eig
  macro plotting(ex)
    return :($(esc(ex)))
  end
  macro isdefined(s::Symbol)
    return isdefined(s)
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

# flatten tuples to arrays
if !@isdefined flatten_tuples
  const flatten_tuples = (x) -> reshape(collect(Iterators.flatten(x)),
                                        length(x[1]), length(x))
end

if !@isdefined ⊗
  const ⊗ = (A,B) -> kron(A, B)
end

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

function transfinite_blend(v1::T, v2, v3, v4, r, s) where T <: Number
  e1 = (α) -> v1 * (1 .- α) / 2 + v3 * (1 .+ α) / 2
  e2 = (α) -> v2 * (1 .- α) / 2 + v4 * (1 .+ α) / 2
  e3 = (α) -> v1 * (1 .- α) / 2 + v2 * (1 .+ α) / 2
  e4 = (α) -> v3 * (1 .- α) / 2 + v4 * (1 .+ α) / 2
  transfinite_blend(e1, e2, e3, e4, r, s)
end
#}}}

#{{{ connectivityarrays
function connectivityarrays(EToV, EToF)
  # number of elements
  nelem = size(EToV, 2)
  nface = maximum(maximum(EToF))

  # Determine secondary arrays
  # FToE : Unique Global Face to Element Number
  # FToLF: Unique Global Face to Element local face number
  # EToO : Element to Unique Global Faces Orientation
  # EToS : Element to Unique Global Face Side

  FToE  = zeros(Int64, 2, nface)
  FToLF = zeros(Int64, 2, nface)
  EToO  = Array{Bool,2}(undef, 4, nelem)
  EToS  = zeros(Int64, 4, nelem)

  # Local Face to Local Vertex map
  LFToLV = flatten_tuples(((1,3), (2, 4), (1,2), (3,4)))
  for e = 1:nelem
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

  er0 = sparse([1  ], [1], [1], Nrp, 1)
  erN = sparse([Nrp], [1], [1], Nrp, 1)
  es0 = sparse([1  ], [1], [1], Nsp, 1)
  esN = sparse([Nsp], [1], [1], Nsp, 1)

  L1 = (Is ⊗ er0')
  L2 = (Is ⊗ erN')
  L3 = (es0' ⊗ Ir)
  L4 = (esN' ⊗ Ir)

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

  τ1 = Diagonal(10./sJ1)
  τ2 = Diagonal(10./sJ2)
  τ3 = Diagonal(10./sJ3)
  τ4 = Diagonal(10./sJ4)


  # TODO: Check signs on Q terms (and update write up with correct signs)
  B1 =  (Sr0 + Sr0') + ((csr0 * Qs + Qs' * csr0) ⊗ Er0) + ((τ1 * H1 * SJ1) ⊗ Er0)
  B2 = -(SrN + SrN') - ((csrN * Qs + Qs' * csrN) ⊗ ErN) + ((τ2 * H2 * SJ2) ⊗ ErN)
  B3 =  (Ss0 + Ss0') + (Es0 ⊗ (crs0 * Qr + Qr' * crs0)) + (Es0 ⊗ (τ3 * H3 * SJ3))
  B4 = -(SsN + SsN') - (EsN ⊗ (crsN * Qr + Qr' * crsN)) + (EsN ⊗ (τ4 * H4 * SJ4))

  F1 =  (Is ⊗ er0') * Sr0 + ((csr0 * Qs) ⊗ er0') + ((τ1 * H1 * SJ1) ⊗ er0')
  F2 = -(Is ⊗ erN') * SrN - ((csrN * Qs) ⊗ erN') + ((τ2 * H2 * SJ2) ⊗ erN')
  F3 =  (es0' ⊗ Ir) * Ss0 + (es0' ⊗ (crs0 * Qr)) + (es0' ⊗ (τ3 * H3 * SJ3))
  F4 = -(esN' ⊗ Ir) * SsN - (esN' ⊗ (crsN * Qr)) + (esN' ⊗ (τ4 * H4 * SJ4))

  @assert B1 ≈ F1' * L1 + L1' * F1 - ((τ1 * H1 * SJ1) ⊗ Er0)
  @assert B2 ≈ F2' * L2 + L2' * F2 - ((τ2 * H2 * SJ2) ⊗ ErN)
  @assert B3 ≈ F3' * L3 + L3' * F3 - (Es0 ⊗ (τ3 * H3 * SJ3))
  @assert B4 ≈ F4' * L4 + L4' * F4 - (EsN ⊗ (τ4 * H4 * SJ4))

  M = A + B1 + B2 + B3 + B4

  # (E, V) = eigen(Matrix(M))
  # println((minimum(E), maximum(E)))
  (M, (F1, F2, F3, F4), (L1, L2, L3, L4), (x, y), Diagonal(J) * (Hs ⊗ Hr),
   (sJ1, sJ2, sJ3, sJ4), (nx1, nx2, nx3, nx4), (ny1, ny2, ny3, ny4),
   (H1I, H2I, H3I, H4I), (τ1, τ2, τ3, τ4))
end
#}}}

let
  #                 1
  #               /   \
  #              1     6
  #             /       \
  #            2    1    3
  #           / \       / \
  #          /   7     9   \
  #         /     \   /     \
  #        /       \ /       \
  #       2         4         5
  #      /    2     |    3     \
  #     /           8           \
  #    /            |            \
  #   5------3------6------4------7
  #
  # verts: Vertices
  verts = ((0,1), (-1/2, 1/2), (1/2, 1/2),
           (0, 1/3), (-1, 0), (0,0), (1,0))

  # EToV: Element to Vertices
  EToV = ((1, 2, 3, 4), (5, 6, 2, 4), (6, 7, 4, 3))

  # EToF: Element to Unique Global Faces
  EToF = ((6, 7, 1, 9), (2, 8, 3, 7), (8, 5, 4, 9))

  # EToN0: Element to sizes
  EToN0 = ((12, 13), (14, 15), (16, 17))

  # FToB: Unique Global Face to Boundary Conditions
  #      -1 = Jumps
  #       0 = internal face
  #       1 = Dirichlet
  #       2 = Neumann
  BC_DIRICHLET        = 1
  BC_NEUMANN          = 2
  BC_LOCKED_INTERFACE = 0
  BC_JUMP_INTERFACE   = -1
  FToB = fill(BC_DIRICHLET, 9)
  FToB[7:9] = BC_NEUMANN

  verts = flatten_tuples(verts)
  EToV  = flatten_tuples(EToV)
  EToF  = flatten_tuples(EToF)
  EToN0 = flatten_tuples(EToN0)

  # number of elements
  nelem = size(EToV, 2)
  nface = size(FToB, 1)

  # Some sanity checks
  @assert typeof(EToV) == Array{Int, 2} && size(EToV) == (4, nelem)
  @assert typeof(EToF) == Array{Int, 2} && size(EToF) == (4, nelem)
  @assert maximum(maximum(EToF)) == nface

  @plotting let
    # Do some plotting
    scatter(verts[1,:], verts[2,:], marker=10, legend=:none)
    for e = 1:nelem
      plot!(verts[1, EToV[[1 2 4 3 1], e]]', verts[2, EToV[[1 2 4 3 1], e]]',
            legend=:none)
    end
    display(plot!())
  end

  # Determine secondary arrays
  # FToE : Unique Global Face to Element Number
  # FToLF: Unique Global Face to Element local face number
  # EToO : Element to Unique Global Faces Orientation
  # EToS : Element to Unique Global Face Side
  (FToE, FToLF, EToO, EToS) = connectivityarrays(EToV, EToF)

  # global mapping
  xg = (r, s)-> r + sin.(π * s) .* cos.(π * r) / 8
  yg = (r, s)-> s - cos.(π * s) .* sin.(π * r) / 8

  # Exact solution
  (kx, ky) = (π, π)
  vex   = (x,y) ->       cos.(kx * x) .* cosh.(ky * y)
  vex_x = (x,y) -> -kx * sin.(kx * x) .* cosh.(ky * y)
  vex_y = (x,y) ->  ky * cos.(kx * x) .* sinh.(ky * y)


  @plotting plot()

  p = 4 # SBP interior order
  ϵ = zeros(3) # size of this array determines the number of levels to run

  OPTYPE = typeof(locoperator(2, 8, 8, (r,s)->r, (r,s)->s))
  for lvl = 1:length(ϵ)
    # println("level = ", lvl)
    lop = Dict{Int64, Tuple{OPTYPE, Array{Float64,1}, UnitRange{Int64}}}()

    Nr = Array{Int64, 1}(undef, nelem)
    Ns = Array{Int64, 1}(undef, nelem)
    Np = Array{Int64, 1}(undef, nelem)
    estart = Array{Int64, 1}(undef, nelem + 1)
    estart[1] = 1
    IM = Array{Int64,1}(undef,0)
    JM = Array{Int64,1}(undef,0)
    VM = Array{Float64,1}(undef,0)
    VH = Array{Float64,1}(undef,0)
    X = Array{Float64,1}(undef,0)
    Y = Array{Float64,1}(undef,0)
    g = Array{Float64,1}(undef,0)
    for e = 1:nelem
      # println("  elm = ", e)

      # Setup number of points in each dimension
      Nr[e] = EToN0[1, e] * (2^(lvl-1))
      Ns[e] = EToN0[2, e] * (2^(lvl-1))
      Np[e] = (Nr[e]+1)*(Ns[e]+1)
      estart[e+1] = estart[e] + Np[e]

      # Setup this elements coordinate transform
      (x1, x2, x3, x4) = verts[1, EToV[:, e]]
      (y1, y2, y3, y4) = verts[2, EToV[:, e]]
      rt = (r,s)->transfinite_blend(x1, x2, x3, x4, r, s)
      st = (r,s)->transfinite_blend(y1, y2, y3, y4, r, s)
      xt = (r,s)->xg(rt(r,s), st(r,s))
      yt = (r,s)->yg(rt(r,s), st(r,s))

      # xt = (r,s)->transfinite_blend(x1, x2, x3, x4, r, s)
      # yt = (r,s)->transfinite_blend(y1, y2, y3, y4, r, s)

      # Build local operators (M assumed Dirichlet BC)
      (M, F, L, (x, y), H, sJ, nx, ny, HfI, τ) =
        locoperator(p, Nr[e], Ns[e], xt, yt)

      @plotting let
        plot!(reshape(x, Nr[e]+1, Ns[e]+1), reshape(y, Nr[e]+1, Ns[e]+1),
              color=:black, legend=:none)
        plot!(reshape(x, Nr[e]+1, Ns[e]+1)', reshape(y, Nr[e]+1, Ns[e]+1)',
              color=:black, legend=:none)
      end

      # Compute the boundary conditions
      v = vex(x,y)
      v_x = vex_x(x,y)
      v_y = vex_y(x,y)

      # Modify operators for the BC
      ge = zeros(Float64, Np[e])
      for f = 1:4
        @plotting plot!(L[f] * x, L[f] * y, color=:red, legend=:none,
                        linewidth=3)
        if FToB[EToF[f, e]] == BC_DIRICHLET
          (xf, yf) = (L[f] * x, L[f] * y)
          vf = vex(xf, yf)
        elseif FToB[EToF[f, e]] == BC_NEUMANN
          (xf, yf) = (L[f] * x, L[f] * y)
          gN = nx[f] .* vex_x(xf, yf) + ny[f] .* vex_y(xf, yf)
          vf = gN ./ diag(τ[f])
          M -= F[f]' * (Diagonal(1 ./ (sJ[f] .* diag(τ[f]))) * HfI[f]) * F[f]
        else
          error("invalid bc")
        end
        ge += F[f]' * vf
      end

      (Ie, Je, Ve) = findnz(M)

      IM = [IM;Ie .+ (estart[e]-1)]
      JM = [JM;Je .+ (estart[e]-1)]
      VM = [VM;Ve]
      g = [g;ge]

      VH = [VH;Vector(diag(H))]
      X = [X;x]
      Y = [Y;y]

      lop[e] = ((M, F, L, (x, y), H, sJ, nx, ny, HfI, τ), ge,
                estart[e]:(estart[e+1]-1))

      #=
      u = M \ ge
      Δ = u - v
      ϵ[lvl] += Δ' * H * Δ
      =#
    end
    # ϵ[lvl] = sqrt(ϵ[lvl])

    @plotting display(plot!())
    TNp = estart[nelem+1]-1
    M = sparse(IM, JM, VM, TNp, TNp)
    H = sparse(1:TNp, 1:TNp, VH, TNp, TNp)
    u = M \ g
    Δ = u - vex(X, Y)
    ϵ[lvl] = sqrt(Δ' * H * Δ)

    println("level = ", lvl, " :: error = ", ϵ[lvl])
  end
  println((log.(ϵ[1:end-1]) - log.(ϵ[2:end])) / log(2))
end

