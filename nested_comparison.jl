include("global_curved.jl")
import PGFPlots
import Metis

Metis.options[Metis.METIS_OPTION_DBGLVL] = Metis.METIS_DBG_SEPINFO

# This bombs out if `N0` and `N1` are set too big
function metis_perm(A)
  OLD_STDOUT = stdout

  Libc.flush_cstdio()
  rdo, wro = redirect_stdout()
  Libc.flush_cstdio()

  perm, _ = Metis.permutation(A)

  Libc.flush_cstdio()
  close(wro)
  Libc.flush_cstdio()
  redirect_stdout(OLD_STDOUT)
  Libc.flush_cstdio()

  lines = readlines(rdo)

  s = zeros(Int64, 4, length(lines))
  for (l, line) in enumerate(lines)
    matches = eachmatch(r"-?\d+\.?\d*", line)
    gen = (parse(Int64, m.match) for m in matches)
    g = collect(gen)
    s[:, l] = g
  end

  (perm, s)
end

function compute_flops(s, n = 1)
  # Flops for the seperator
  flops = s[4, n]^3
  if size(s, 2) > n && s[2, n] == s[1, n + 1]
    # @show "branch"
    # Flops for the right branch
    (n, flops_right) = compute_flops(s, n + 1)
    # Flops for the left branch
    (n, flops_left ) = compute_flops(s, n + 1)
    flops += flops_right + flops_left
  else
    # @show "leaf"
    # For the leaf case the two sides are just the flops we need
    flops += s[2, n]^3 + s[3, n]^3
  end
  # @show (n, flops)
  (n, flops)
end

let
  # SBP interior order
  SBPp   = 6

  # mesh file side set type to actually boundary condition type
  bc_map = [BC_DIRICHLET, BC_DIRICHLET, BC_NEUMANN, BC_NEUMANN,
            BC_JUMP_INTERFACE]
  (verts, EToV, EToF, FToB, EToDomain) = read_inp_2d("meshes/square_circle.inp";
                                                     bc_map = bc_map)
  # EToV defines the element by its vertices
  # EToF defines element by its four faces, in global face number
  # FToB defines whether face is Dirichlet (1), Neumann (2), interior jump (7)
  #      or just an interior interface (0)
  # EToDomain is 1 if element is inside circle; 2 otherwise

  # number of elements and faces
  (nelems, nfaces) = (size(EToV, 2), size(FToB, 1))
  @show (nelems, nfaces)

  # This is needed to fix up points that should be on the boundary of the
  # circle, but the mesh didn't quite put them there
  for v in 1:size(verts, 2)
    x, y = verts[1, v], verts[2, v]
    if abs(hypot(x,y) - 1) < 1e-5
      Q = atan(y, x)
      verts[1, v], verts[2, v] = cos(Q), sin(Q)
    end
  end

  # Plot the original connectivity before mesh warping
  # plot_connectivity(verts, EToV)

  # Increasing this size is how to get it to break...
  N1 = N0 =  17

  # EToN0 is the base mesh size (e.g., before refinement)
  EToN0 = zeros(Int64, 2, nelems)
  EToN0[1, :] .= N0
  EToN0[2, :] .= N1

  @assert typeof(EToV) == Array{Int, 2} && size(EToV) == (4, nelems)
  @assert typeof(EToF) == Array{Int, 2} && size(EToF) == (4, nelems)
  @assert maximum(maximum(EToF)) == nfaces

  # Determine secondary arrays
  # FToE : Unique Global Face to Element Number
  #        (the i'th column of this stores the element numbers that share the
  #        global face number i)
  # FToLF: Unique Global Face to Element local face number
  #        (the i'th column of this stores the element local face numbers that
  #        shares the global face number i)
  # EToO : Element to Unique Global Faces Orientation
  #        (the i'th column of this stores the whether the element and global
  #        face are oriented in the same way in physical memory or need to be
  #        rotated)
  # EToS : Element to Unique Global Face Side
  #        (the i'th column of this stores whether an element face is on the
  #        plus side or minus side of the global face)
  (FToE, FToLF, EToO, EToS) = connectivityarrays(EToV, EToF)

  # Exact solution
  Lx = maximum(verts[1,:])
  Ly = maximum(abs.(verts[2,:]))

  # Set up the local grid dimensions
  Nr = EToN0[1, :]
  Ns = EToN0[2, :]

  #
  # Build the local volume operators
  #

  # Dictionary to store the operators
  OPTYPE = typeof(locoperator(2, 16, 16))
  lop = Dict{Int64, OPTYPE}()

  # Loop over blocks and create local operators
  for e = 1:nelems
    # Get the element corners
    (x1, x2, x3, x4) = verts[1, EToV[:, e]]
    (y1, y2, y3, y4) = verts[2, EToV[:, e]]

    # Initialize the block transformations as transfinite between the corners
    ex = [(α) -> x1 * (1 .- α) / 2 + x3 * (1 .+ α) / 2,
          (α) -> x2 * (1 .- α) / 2 + x4 * (1 .+ α) / 2,
          (α) -> x1 * (1 .- α) / 2 + x2 * (1 .+ α) / 2,
          (α) -> x3 * (1 .- α) / 2 + x4 * (1 .+ α) / 2]
    exα = [(α) -> -x1 / 2 + x3 / 2,
           (α) -> -x2 / 2 + x4 / 2,
           (α) -> -x1 / 2 + x2 / 2,
           (α) -> -x3 / 2 + x4 / 2]
    ey = [(α) -> y1 * (1 .- α) / 2 + y3 * (1 .+ α) / 2,
          (α) -> y2 * (1 .- α) / 2 + y4 * (1 .+ α) / 2,
          (α) -> y1 * (1 .- α) / 2 + y2 * (1 .+ α) / 2,
          (α) -> y3 * (1 .- α) / 2 + y4 * (1 .+ α) / 2]
    eyα = [(α) -> -y1 / 2 + y3 / 2,
           (α) -> -y2 / 2 + y4 / 2,
           (α) -> -y1 / 2 + y2 / 2,
           (α) -> -y3 / 2 + y4 / 2]

    # For blocks on the circle, put in the curved edge transform
    if FToB[EToF[1, e]] == BC_JUMP_INTERFACE
      error("curved face 1 not implemented yet")
    end
    if FToB[EToF[2, e]] == BC_JUMP_INTERFACE
      error("curved face 2 not implemented yet")
    end
    if FToB[EToF[3, e]] == BC_JUMP_INTERFACE
      Q1 = atan(y1, x1)
      Q2 = atan(y2, x2)
      if !(-π/2 < Q1 - Q2 < π/2)
        Q2 -= sign(Q2) * 2 * π
      end
      ex[3] = (α) -> cos.(Q1 * (1 .- α) / 2 + Q2 * (1 .+ α) / 2)
      ey[3] = (α) -> sin.(Q1 * (1 .- α) / 2 + Q2 * (1 .+ α) / 2)
      β3 = (Q2 - Q1) / 2
      exα[3] = (α) -> -β3 .* sin.(Q1 * (1 .- α) / 2 + Q2 * (1 .+ α) / 2)
      eyα[3] = (α) -> +β3 .* cos.(Q1 * (1 .- α) / 2 + Q2 * (1 .+ α) / 2)
    end
    if FToB[EToF[4, e]] == BC_JUMP_INTERFACE
      Q3 = atan(y3, x3)
      Q4 = atan(y4, x4)
      if !(-π/2 < Q3 - Q4 < π/2)
        error("curved face 4 angle correction not implemented yet")
      end
      ex[4] = (α) -> cos.(Q3 * (1 .- α) / 2 + Q4 * (1 .+ α) / 2)
      ey[4] = (α) -> sin.(Q3 * (1 .- α) / 2 + Q4 * (1 .+ α) / 2)
      β4 = (Q4 - Q3) / 2
      exα[4] = (α) -> -β4 .* sin.(Q3 * (1 .- α) / 2 + Q4 * (1 .+ α) / 2)
      eyα[4] = (α) -> +β4 .* cos.(Q3 * (1 .- α) / 2 + Q4 * (1 .+ α) / 2)
    end

    # Create the volume transform as the transfinite blending of the edge
    # transformations
    xt(r,s) = transfinite_blend(ex[1], ex[2], ex[3], ex[4],
                                exα[1], exα[2], exα[3], exα[4],
                                r, s)
    yt(r,s) = transfinite_blend(ey[1], ey[2], ey[3], ey[4],
                                eyα[1], eyα[2], eyα[3], eyα[4],
                                r, s)

    metrics = create_metrics(SBPp, Nr[e], Ns[e], xt, yt)

    # Build local operators
    lop[e] = locoperator(SBPp, Nr[e], Ns[e], metrics, FToB[EToF[:, e]])
  end

  # If this is the first mesh level plot the mesh
  plot_blocks(lop)

  #
  # Do some assemble of the global volume operators
  #
  (M, FbarT, D, vstarts, FToλstarts) =
  LocalGlobalOperators(lop, Nr, Ns, FToB, FToE, FToLF, EToO, EToS,
                       (x) -> cholesky(Symmetric(x)))
  locfactors = M.F

  # Get a unique array indexes for the face to jumps map
  FToδstarts = bcstarts(FToB, FToE, FToLF, BC_JUMP_INTERFACE, Nr, Ns)

  # Compute the number of volume, trace (λ), and jump (δ) points
  VNp = vstarts[nelems+1]-1
  λNp = FToλstarts[nfaces+1]-1
  δNp = FToδstarts[nfaces+1]-1
  # @show (VNp, λNp, δNp)

  # trace system
  B = assembleλmatrix(FToλstarts, vstarts, EToF, FToB, locfactors, D, FbarT)
  B2 = (B + B')/2
  @assert B2 ≈ B
  (perm_B, s_B) = metis_perm(B2)
  B2_p = B2[perm_B, perm_B]

  # Monolithic system
  M = blockdiag(ntuple(i->lop[i].M̃, length(lop))...)
  A = [M FbarT'; FbarT Diagonal(D)]
  A2 = (A + A')/2
  @assert A2 ≈ A
  (perm_A, s_A) = metis_perm(A2)
  A2_p = A2[perm_A, perm_A]

  # Displacement system
  M = blockdiag(ntuple(i->lop[i].M̃, length(lop))...)
  C = M - FbarT' * Diagonal(1 ./ D) * FbarT
  C2 = (C + C')/2
  @assert C2 ≈ C
  (perm_C, s_C) = metis_perm(C2)
  C2_p = C2[perm_C, perm_C]

  # Single block Displacement system
  Np = (N1 + 1)*(N0 + 1)
  G = M[1:Np, 1:Np]
  G2 = (G + G')/2
  @assert G2 ≈ G
  (perm_G, s_G) = metis_perm(G2)
  G2_p = G2[perm_G, perm_G]

  println("flops for trace system inverse:       $(compute_flops(s_B)[2])")
  println("flops for monolithic system inverse:  $(compute_flops(s_A)[2])")
  println("flops for displacement system inverse $(compute_flops(s_C)[2])")
  println("flops for single block inverse:       $(compute_flops(s_G)[2])")
  println("flops for all single block inverse:   $(nelems * compute_flops(s_G)[2])")
  println("flops for hybridized inverse:         $(nelems * compute_flops(s_G)[2] + compute_flops(s_B)[2])")
end
nothing
