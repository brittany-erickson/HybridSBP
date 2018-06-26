include("diagonal_sbp.jl")
using SparseArrays

function make_operators_dirchlet(p, Nx, Ny)
  (D2x_1d, BSx_1d, HIx_1d, Hx_1d, rx_1d) = diagonal_sbp_D2(p, Nx)
  (D2y_1d, BSy_1d, HIy_1d, Hy_1d, ry_1d) = diagonal_sbp_D2(p, Ny)

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

  Bn = I - In * In'

  x = kron(ones(Int64, Ny+1), rx_1d)
  y = kron(ry_1d, ones(Int64, Nx+1))

  (Ax, Ay, H, In, Bn, x, y)
end


#{{{ Build Macros element maps for continuous corners and faces
function buildmacromaps(EToE, EToF, EToO)
  @assert size(EToE) == size(EToF) == size(EToO)
  K = size(EToE, 2) # number of macro elements

  #{{{ Element to Continuous Face
  nfaces = 4
  EToCF = zeros(Int64, nfaces, K)
  gf = 0
  for e = 1:K
    for f = 1:nfaces
      if EToE[f, e] == 0
        gf += 1
        EToCF[f, e] = gf
      elseif EToE[f, e] < e
        ne = EToE[f, e]
        nf = EToF[f, e]
        EToCF[f, e] = EToCF[nf, ne]
      else
        gf += 1
        EToCF[f, e] = gf
      end
    end
  end
  CFToE = zeros(Int64, 2, gf)
  CFToF = zeros(Int64, 2, gf)
  CFToO = zeros(Int64, 2, gf)
  for e = 1:K
    for f = 1:nfaces
      gf = EToCF[f, e]
      if CFToE[1, gf] == 0
        CFToE[1, gf] = e
        CFToF[1, gf] = f
        CFToO[2, gf] = 0
      else
        @assert CFToE[2, gf] == 0
        CFToE[2, gf] = e
        CFToF[2, gf] = f
        CFToO[2, gf] = EToO[f, e]
      end
    end
  end
  #}}}

  #{{{ Element to Continuous Face
  ncorners = 4
  EToECE = zeros(Int64, ncorners, K) # element to element who owns continuous elem
  EToECC = zeros(Int64, ncorners, K) # element to element who owns continuous corn

  # First pass
  FToC = [1 2 1 3; 3 4 2 4]; 
  for e = 1:K
    (e1, e2, e3, e4) = EToE[:, e]
    e1 = (e1 > 0) ? e1 : e
    e2 = (e2 > 0) ? e2 : e
    e3 = (e3 > 0) ? e3 : e
    e4 = (e4 > 0) ? e4 : e

    (f1, f2, f3, f4) = EToF[:, e]
    (o1, o2, o3, o4) = EToO[:, e]

    # c1 -> (f1, f3)
    c = 1
    EToECE[c, e] = em = min(e, e1, e3)
    if e == em
      EToECC[c, e] = c
    elseif e1 < e3
      EToECC[c, e] = (o1 == 0) ? FToC[1, f1] : -FToC[2, f1]
    else
      EToECC[c, e] = (o3 == 0) ? FToC[1, f3] : FToC[2, f3]
    end

    # c2 -> (f2, f3)
    c = 2
    EToECE[c, e] = em = min(e, e2, e3)
    if e == em
      EToECC[c, e] = c
    elseif e2 < e3
      EToECC[c, e] = (o2 == 0) ? FToC[1, f2] : FToC[2, f2]
    else
      EToECC[c, e] = (o3 == 0) ? FToC[2, f3] : FToC[1, f3]
    end

    # c3 -> (f1, f4)
    c = 3
    EToECE[c, e] = em = min(e, e1, e4)
    if e == em
      EToECC[c, e] = c
    elseif e1 < e4
      EToECC[c, e] = (o1 == 0) ? FToC[2, f1] : FToC[1, f1]
    else
      EToECC[c, e] = (o4 == 0) ? FToC[1, f4] : FToC[2, f4]
    end

    # c4 -> (f2, f4)
    c = 4
    EToECE[c, e] = em = min(e, e2, e4)
    if e == em
      EToECC[c, e] = c
    elseif e2 < e4
      EToECC[c, e] = (o2 == 0) ? FToC[2, f2] : FToC[1, f2]
    else
      EToECC[c, e] = (o4 == 0) ? FToC[2, f4] : FToC[1, f4]
    end
  end

  changes = 1
  while changes > 0
    changes = 0
    for e = 1:K
      for c = 1:ncorners
        ne = EToECE[c, e]
        nc = EToECC[c, e]
        nne = EToECE[nc, ne]
        nnc = EToECC[nc, ne]
        # If my neighbors neighbor is different than me, I might need to update
        if (ne, nc) != (nne, nnc)
          # If my neighbor elm number is smaller, use neighbor value
          if nne < ne
            (EToECE[c, e], EToECC[c, e]) = (nne, nnc)
            changes += 1
          # If my neighbor elm number is mine, but corner number is lower use
          # that one
          else nne == ne && nnc < nc
            EToECC[c, e] = nnc
            changes += 1
          end
        end
      end
    end
  end

  # Build element corner to continous corner map
  EToCC = zeros(Int64, ncorners, K)
  nglobalcorners = 0
  for e = 1:K
    for c = 1:ncorners
      ne = EToECE[c, e]
      nc = EToECC[c, e]
      if ne == e && nc == c
        nglobalcorners += 1
        EToCC[c, e] = nglobalcorners
      else
        EToCC[c, e] = EToCC[nc, ne]
      end
    end
  end

  # Build continuous corner to element maps
  # Continuous corner multiplicity
  CCMultiplicity = zeros(Int64, nglobalcorners)
  for e = 1:K
    for c = 1:ncorners
      CCMultiplicity[EToCC[c, e]] += 1
    end
  end
  # Continuous corner offsets
  CCOffsets = [0;accumulate(+, CCMultiplicity)]

  # Continuous corner To element map
  CCToE = zeros(Int64, K * ncorners)

  # Continuous corner To element corner map
  CCToC = zeros(Int64, K * ncorners)
  for e = 1:K
    for c = 1:ncorners
      nc = EToCC[c, e]
      offset = CCOffsets[nc]+1
      nextoffset = CCOffsets[nc+1]+1
      while CCToE[offset] != 0
        offset += 1
        @assert offset < nextoffset
      end
      CCToE[offset] = e
      CCToC[offset] = c
    end
  end
  #}}}

  (EToCC, EToCF, CFToE, CFToF, CFToO, CCOffsets, CCToE, CCToC)
end
#}}}

#{{{ Build Macros element maps for continuous corners and faces
function buildCNToDN(EToN, EToCC, EToCF, CFToE, CFToF, CFToO)
  K = size(EToN, 2)
  num_unique_corners = maximum(EToCC)
  num_unique_faces   = size(CFToE, 2)

  face_offsets = zeros(Int64, num_unique_faces+1)
  face_offsets[1] = num_unique_corners+1
  for gf = 1:num_unique_faces
    e1 = CFToE[1, gf]

    # Number of points along this face
    Nf = EToN[2-div(CFToF[1,gf]-1,2), e1]

    # make sure boundary or sizes match along face
    e2 = CFToE[2, gf]
    @assert e2 == 0 || Nf == EToN[2-div(CFToF[2,gf]-1,2), e2]

    face_offsets[gf+1] = face_offsets[gf] + (Nf - 1)
  end

  num_corn_points = num_unique_corners
  num_face_points = face_offsets[end] - face_offsets[1]
  num_volm_points = sum(prod(EToN .- 1, dims=1))
  num_cont_points = num_corn_points + num_face_points + num_volm_points
  num_disc_points = sum(prod(EToN .+ 1, dims=1))

  # IC2D: dontinuous to disconcinous
  InI = 1:num_disc_points
  D2C = zeros(Int64, num_disc_points)
  InV = ones(Int64, num_disc_points)
  vs = num_corn_points + num_face_points
  ds = 0
  for e = 1:K
    (gc1, gc2, gc3, gc4) = EToCC[:, e]
    (gf1, gf2, gf3, gf4) = EToCF[:, e]
    o1 = CFToO[(CFToE[1,gf1] == e && CFToF[1,gf1] == 1) ? 1 : 2, gf1]
    o2 = CFToO[(CFToE[1,gf2] == e && CFToF[1,gf2] == 1) ? 1 : 2, gf2]
    o3 = CFToO[(CFToE[1,gf3] == e && CFToF[1,gf3] == 1) ? 1 : 2, gf3]
    o4 = CFToO[(CFToE[1,gf4] == e && CFToF[1,gf4] == 1) ? 1 : 2, gf4]

    fp1 = (o1 == 0) ? (face_offsets[gf1]:face_offsets[gf1+1]-1) :
                      (face_offsets[gf1+1]-1:-1:face_offsets[gf1])
    fp2 = (o2 == 0) ? (face_offsets[gf2]:face_offsets[gf2+1]-1) :
                      (face_offsets[gf2+1]-1:-1:face_offsets[gf2])
    fp3 = (o3 == 0) ? (face_offsets[gf3]:face_offsets[gf3+1]-1) :
                      (face_offsets[gf3+1]-1:-1:face_offsets[gf3])
    fp4 = (o4 == 0) ? (face_offsets[gf4]:face_offsets[gf4+1]-1) :
                      (face_offsets[gf4+1]-1:-1:face_offsets[gf4])

    #
    vp = vs .+ (1:(EToN[1,e]-1) * (EToN[2,e]-1))
    dp = ds .+ (1:(EToN[1,e]+1) * (EToN[2,e]+1))

    D2C[dp] = ([ gc1 collect(fp1)' gc3;
             collect(fp3) reshape(vp, EToN[1,e]-1, EToN[2,e]-1) collect(fp4);
             gc2 collect(fp2)' gc4])

    vs = vp[end]
    ds = dp[end]
  end
  IC2D = sparse(InI, D2C, InV)

  #IB2C: Boundary to Continuous
  #II2C: Interior to Continuous
  num_cf_points = num_corn_points+num_face_points
  Bpts = zeros(Int64, num_cf_points)
  Doffsets = [0; accumulate(+, prod(EToN .+ 1; dims=1)[:])]
  for gf = 1:num_unique_faces
    # If this face has only 1 side, then we're a boundary face
    if CFToE[2, gf] == 0
      e = CFToE[1, gf]
      f = CFToF[1, gf]
      offset = Doffsets[e]
      (Npx, Npy) = EToN[:, e] .+ 1
      if f == 1
        Bpts[D2C[offset .+ (1:Npx:(Npx*Npy))]] .= 1
      elseif f == 2
        Bpts[D2C[offset .+ (Npx:Npx:(Npx*Npy))]] .= 1
      elseif f == 3
        Bpts[D2C[offset .+ (1:Npx)]] .= 1
      elseif f == 4
        Bpts[D2C[offset .+ Npx * (Npy-1) .+ (1:Npx)]] .= 1
      end
    end
  end
  nbcs = sum(Bpts)

  B2C = zeros(Int64, nbcs)
  I2C = zeros(Int64, num_cont_points-nbcs)
  b = 0
  i = 0
  for n = 1:num_cf_points
    if Bpts[n] == 1
      b += 1
      B2C[b] = n
    else
      i += 1
      I2C[i] = n
    end
  end
  I2C[i+1:num_cont_points-nbcs] = (num_cf_points+1):num_cont_points
  II2C = sparse(I2C, 1:num_cont_points-nbcs, ones(Int64, num_cont_points-nbcs),
                num_cont_points, num_cont_points-nbcs)
  IB2C = sparse(B2C, 1:nbcs, ones(Int64, nbcs), num_cont_points, nbcs)

  (IC2D, IB2C, II2C)
end
#}}}

#
# +---+---+
# | 3 | 4 |
# +---+---+
# | 1 | 2 |
# +---+---+
let
  EToE = [0 2 0 3;
          1 0 0 4;
          0 4 1 0;
          3 0 2 0]';
  EToF = [0 1 0 3;
          2 0 0 3;
          0 1 4 0;
          2 0 4 0]';
  EToO = [0 0 0 0;
          0 0 0 0;
          0 0 0 0;
          0 0 0 0]';

  ϵ = zeros(3)
  N0 = (2, 8, 12, 16, 21)
  for p = 2:2:10
    for j = 0:length(ϵ)-1
      Nx1 = 2^j * (N0[div(p,2)] + 0)
      Nx2 = 2^j * (N0[div(p,2)] + 1)
      Ny1 = 2^j * (N0[div(p,2)] + 2)
      Ny2 = 2^j * (N0[div(p,2)] + 3)


      EToN = [Nx1 Nx2 Nx1 Nx2;
              Ny1 Ny1 Ny2 Ny2]

      (EToCC, EToCF, CFToE, CFToF, CFToO, ~, ~, ~) = buildmacromaps(EToE, EToF, EToO)

      (IC2D, IB2C, II2C) = buildCNToDN(EToN, EToCC, EToCF, CFToE, CFToF, CFToO)
      D2C_scaling = 1 ./ sum(IC2D, dims=1)[:] # scaling for going from D to C

      (Ax1, Ay1, H1, ~,~, x1, y1) = make_operators_dirchlet(p, EToN[1,1], EToN[2,1])
      x1 .-= 1
      y1 .-= 1
      (Ax2, Ay2, H2, ~,~, x2, y2) = make_operators_dirchlet(p, EToN[1,2], EToN[2,2])
      x2 .+= 1
      y2 .-= 1
      (Ax3, Ay3, H3, ~,~, x3, y3) = make_operators_dirchlet(p, EToN[1,3], EToN[2,3])
      x3 .-= 1
      y3 .+= 1
      (Ax4, Ay4, H4, ~,~, x4, y4) = make_operators_dirchlet(p, EToN[1,4], EToN[2,4])
      x4 .+= 1
      y4 .+= 1

      # Continuous x and y
      xc = D2C_scaling .* (IC2D' * [x1;x2;x3;x4])
      yc = D2C_scaling .* (IC2D' * [y1;y2;y3;y4])
      @assert IC2D * xc ≈ [x1;x2;x3;x4]
      @assert IC2D * yc ≈ [y1;y2;y3;y4]

      A1 = Ax1 + Ay1
      A2 = Ax2 + Ay2
      A3 = Ax3 + Ay3
      A4 = Ax4 + Ay4

      A = cat(A1, A2, A3, A4, dims=[1,2])
      H = cat(H1, H2, H3, H4, dims=[1,2])
      AC = IC2D' * A * IC2D
      ACI = II2C' * AC * II2C

      ue = sin.(π * xc) .* sinh.(π * yc)
      ub = IB2C' * ue
      uhi = ACI \ (-II2C' * (AC * (IB2C * ub)))
      uh = (IB2C * ub) + (II2C * uhi)

      Δu = IC2D * (uh - ue)
      ϵ[j+1] = sqrt(Δu' * H * Δu)
      println((p, j, ϵ[j+1]))
    end
    println((log.(ϵ[1:end-1]) - log.(ϵ[2:end])) / log(2))
    println()
  end
end
