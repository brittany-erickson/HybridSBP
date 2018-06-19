include("diagonal_sbp_D2.jl")
using SparseArrays

function locoperator(p, Nx, Ny, corners::NTuple{4, Tuple{T, T}}) where T <: Number
  @assert (corners[1][1], corners[2][1]) == (corners[3][1], corners[4][1])
  @assert (corners[1][2], corners[3][2]) == (corners[2][2], corners[4][2])

  (D2x_1d, BSx_1d, HIx_1d, Hx_1d, rx_1d) = diagonal_sbp_D2(p, Nx; xc=
                                                           (corners[1][1],
                                                            corners[4][1]))
  (D2y_1d, BSy_1d, HIy_1d, Hy_1d, ry_1d) = diagonal_sbp_D2(p, Ny; xc=
                                                           (corners[1][2],
                                                            corners[4][2]))

  Ax_1d = BSx_1d - Hx_1d * D2x_1d
  Ay_1d = BSy_1d - Hy_1d * D2y_1d

  Ax = kron(Hy_1d, Ax_1d)
  Ay = kron(Ay_1d, Hx_1d)

  H  = kron(Hy_1d, Hx_1d)

  E0x_1d = sparse([     1], [     1], [1.], Nx + 1, Nx + 1)
  ENx_1d = sparse([Nx + 1], [Nx + 1], [1.], Nx + 1, Nx + 1)
  E0y_1d = sparse([     1], [     1], [1.], Ny + 1, Ny + 1)
  ENy_1d = sparse([Ny + 1], [Ny + 1], [1.], Ny + 1, Ny + 1)

  Inx = sparse(2:Nx,1:Nx-1, ones(Int64, Nx-1), Nx+1, Nx-1)
  Iny = sparse(2:Ny,1:Ny-1, ones(Int64, Ny-1), Ny+1, Ny-1)
  In = kron(Iny, Inx)

  Nfp = 2 * Nx + 2 * Ny
  Nxp = Nx + 1
  Nyp = Ny + 1
  B2C = [1; Nxp; 1+Nxp*Ny; Nxp*Nyp;
         ((Nxp+1):Nxp:Nxp*(Ny-1)+1);
         ((Nxp+1):Nxp:Nxp*(Ny-1)+1) .+ Nx;
         (2:Nx);
         (2:Nx) .+ Nxp*Ny]
  Bn = sparse(B2C, 1:Nfp, ones(Int64, Nfp), Nxp * Nyp, Nfp)

  x = kron(ones(Int64, Ny+1), rx_1d)
  y = kron(ry_1d, ones(Int64, Nx+1))

  (Ax, Ay, H, In, Bn, x, y)
end

function locsolve(Ax, Ay, In, Bn, uB)
  A = Ax + Ay
  AI = (In' * A * In)
  AB = In' * A * Bn
  uI = AI \ (-AB * uB)
end

let
  p = 8
  N0 = (2, 8, 12, 16, 21)
  ϵ = zeros(4)
  for lvl = 1:length(ϵ)
    Nx1 = 2^(lvl-1) * (N0[div(p,2)] + 0)
    Nx2 = 2^(lvl-1) * (N0[div(p,2)] + 1)
    Ny1 = 2^(lvl-1) * (N0[div(p,2)] + 2)
    Ny2 = 2^(lvl-1) * (N0[div(p,2)] + 3)

    # 7---9---8--13---9
    # |       |       |
    # 4   3   5   4   6
    # |       |       |
    # 4---8---5--12---6
    # |       |       |
    # 1   1   2   2   3
    # |       |       |
    # 1---7---2--11---3
    # v: Vertices
    v = ((-1,-1), ( 0,-1), ( 1,-1),
         (-1, 0), ( 0, 0), ( 1, 0),
         (-1, 1), ( 0, 1), ( 1, 1))

    # EToV: Element to Vertices
    EToV = ((1, 2, 4, 5), (2, 3,  5,  6), (4, 5, 7, 8), (5, 6,  8,  9))

    # EToF: Element to Unique Global Faces
    EToF = ((1, 2, 7, 8), (2, 3, 10, 11), (4, 5, 8, 9), (5, 6, 11, 12))

    # EToF: Element to Unique Global Faces Orientation
    EToO = ((0, 0, 0, 0), (0, 0,  0,  0), (0, 0, 0, 0), (0, 0,  0,  0))

    # FToB: Unique Global Face to Boundary Conditions
    #       0 = internal face
    #       1 = Dirichlet
    #       2 = Neumann (not supported yet)
    FToB = (1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1)

    # EToN: Element to (Nx, Ny) size
    #       (N[xy[ = number of points in that dimension minus 1)
    EToN = ((Nx1, Ny1), (Nx2, Ny1), (Nx1, Ny2), (Nx2, Ny2))

    # FToNp: Face to Face Number of points (does not include the corners)
    FToNp = zeros(Int64, length(FToB))
    for e = 1:length(EToV)
      for lf = 1:length(EToV[e])
        gf = EToF[e][lf]
        if FToNp[gf] == 0
          FToNp[gf] = EToN[e][2 - div(lf - 1, 2)] - 1
        else
          @assert FToNp[gf] == EToN[e][2 - div(lf - 1, 2)] - 1
        end
      end
    end

    # FToOffset: Face to Face Offset in unique global numbering of trace
    FToOffset = accumulate(+, [length(v) + 1; FToNp[:]])

    num_elm = length(EToV)
    Ax = Array{SparseMatrixCSC{Float64, Int64}, 1}(undef, num_elm)
    Ay = Array{SparseMatrixCSC{Float64, Int64}, 1}(undef, num_elm)
    H  = Array{SparseMatrixCSC{Float64, Int64}, 1}(undef, num_elm)
    In = Array{SparseMatrixCSC{  Int64, Int64}, 1}(undef, num_elm)
    Bn = Array{SparseMatrixCSC{  Int64, Int64}, 1}(undef, num_elm)
    x  = Array{Array{Float64, 1}, 1}(undef, num_elm)
    y  = Array{Array{Float64, 1}, 1}(undef, num_elm)

    # unique point [xy]trace along the elements
    xtrace = fill(NaN, FToOffset[end]-1)
    ytrace = fill(NaN, FToOffset[end]-1)
    for e = 1:num_elm
      (Ax[e], Ay[e], H[e], In[e], Bn[e], x[e], y[e]) =
        locoperator(p, EToN[e][1], EToN[e][2], (v[EToV[e][1]], v[EToV[e][2]],
                                                v[EToV[e][3]], v[EToV[e][4]]))
      xet = Bn[e]' * x[e]
      yet = Bn[e]' * y[e]
      for c = 1:length(EToV[e])
        if isnan(xtrace[EToV[e][c]]) && isnan(ytrace[EToV[e][c]])
          (xtrace[EToV[e][c]], ytrace[EToV[e][c]]) = (xet[c], yet[c])
        else
          @assert xtrace[EToV[e][c]] ≈ xet[c]
          @assert ytrace[EToV[e][c]] ≈ yet[c]
        end
      end

      st = length(EToV[e])
      for lf = 1:length(EToF[e])
        gf = EToF[e][lf]
        glorng = FToOffset[gf]:(FToOffset[gf+1]-1)
        locrng = st .+ (1:length(glorng))
        if isnan(xtrace[glorng[1]])
          xtrace[glorng] = xet[locrng]
          ytrace[glorng] = yet[locrng]
        else
          @assert xtrace[glorng] ≈ xet[locrng]
          @assert ytrace[glorng] ≈ yet[locrng]
        end
        st = locrng[end]
      end
    end
    @assert !maximum(isnan.(xtrace)) && !maximum(isnan.(ytrace))

    uexact = (x,y) -> sin.(π * x) .* sinh.(π * y)
    utrace = uexact(xtrace, ytrace)
    ϵ[lvl] = 0
    for e = 1:num_elm
      (c1, c2, c3, c4) = EToV[e]
      (gf1, gf2, gf3, gf4) = EToF[e]
      f1 = FToOffset[gf1]:(FToOffset[gf1+1]-1)
      f2 = FToOffset[gf2]:(FToOffset[gf2+1]-1)
      f3 = FToOffset[gf3]:(FToOffset[gf3+1]-1)
      f4 = FToOffset[gf4]:(FToOffset[gf4+1]-1)
      uB = utrace[[c1;c2;c3;c4;f1;f2;f3;f4]]
      uI = locsolve(Ax[e], Ay[e], In[e], Bn[e], uB)
      uh = In[e] * uI + Bn[e] * uB
      Δu = uh - uexact(x[e], y[e])
      ϵ[lvl] += Δu' * H[e] * Δu
    end
    ϵ[lvl] = sqrt(ϵ[lvl])
    println("level = ", lvl, " error = ", ϵ[lvl])
  end
  println((log.(ϵ[1:end-1]) - log.(ϵ[2:end])) / log(2))
  println()

  nothing
end
