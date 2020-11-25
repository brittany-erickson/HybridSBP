include("global_curved.jl")
import PGFPlots

let
  method = :displacements

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

  # This is the base mesh size in each dimension
  N1 = N0 = 17

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

  pgf_axis = PGFPlots.Axis(style="width=5cm, height=5cm, ticks=none",
                           xlabel=PGFPlots.L"$x$",
                           ylabel=PGFPlots.L"$y$",
                           xmin = -2, xmax = 2,
                           ymin = -2, ymax = 2)
  for f in 1:nfaces
    if FToB[f] != BC_JUMP_INTERFACE
      (e, lf) = FToE[1,f], FToLF[1,f]
      if lf == 1
        v1, v2 = EToV[1, e], EToV[3, e]
      elseif lf == 2
        v1, v2 = EToV[2, e], EToV[4, e]
      elseif lf == 3
        v1, v2 = EToV[1, e], EToV[2, e]
      else
        v1, v2 = EToV[3, e], EToV[4, e]
      end
      x = verts[1, [v1 v2]][:]
      y = verts[2, [v1 v2]][:]
      push!(pgf_axis, PGFPlots.Linear(x, y, style="no marks, solid, black"))
    end
  end
  push!(pgf_axis, PGFPlots.Circle(0,0,1, style = "very thick, red"))
  PGFPlots.save("square_circle.tikz", pgf_axis)

  # Exact solution
  Lx = maximum(verts[1,:])
  Ly = maximum(abs.(verts[2,:]))



  c = exp(1)/(1 + exp(1))

  vinside(x,y,e) = begin
    r = sqrt.(x .^ 2 + y .^ 2)
    theta = atan.(y,x)
    return c .* (1 .- exp.(-1 .* r .^ 2)) .* r .* sin.(theta)
  end

  voutside(x,y,e) = begin
    r = sqrt.(x .^ 2 + y .^ 2)
    theta = atan.(y,x)
    return (r .- 1) .^ 2 .* cos.(theta) .+ (r .- 1) .* sin.(theta)
  end


  vinside_x(x,y,e) = begin
    r = sqrt.(x .^ 2 + y .^ 2)
    theta = atan.(y,x)
    dtheta_dx = -1 .* sin.(theta) ./ r
    dr_dx = cos.(theta)
    dv_dr = c * (2 .* r .^ 2 .* exp.(-1 .* r .^ 2) .+ 1 .- exp.(-1 .* r .^ 2)) .* sin.(theta)
    dv_dtheta = c .* (1 .- exp.(-1 .* r .^ 2)) .* r .* cos.(theta)
    return dv_dr .* dr_dx + dv_dtheta .* dtheta_dx
  end


  vinside_y(x,y,e) = begin
    r = sqrt.(x .^ 2 + y .^ 2)
    theta = atan.(y,x)
    dtheta_dy = cos.(theta) ./ r
    dr_dy = sin.(theta)
    dv_dr = c * (2 .* r .^ 2 .* exp.(-1 .* r .^ 2) .+ 1 .- exp.(-1 .* r .^ 2)) .* sin.(theta)
    dv_dtheta = c .* (1 .- exp.(-1 .* r .^ 2)) .* r .* cos.(theta)
    return dv_dr .* dr_dy + dv_dtheta .* dtheta_dy
  end

  voutside_x(x,y,e) = begin
    r = sqrt.(x .^ 2 + y .^ 2)
    theta = atan.(y,x)
    dtheta_dx = -1 .* sin.(theta) ./ r
    dr_dx = cos.(theta)
    dv_dr =  2 .* (r .- 1) .* cos.(theta) .+ sin.(theta)
    dv_dtheta = -1 .* (r .- 1) .^ 2 .* sin.(theta) .+ (r .- 1) .* cos.(theta)
    return dv_dr .* dr_dx + dv_dtheta .* dtheta_dx
  end


  voutside_y(x,y,e) = begin
    r = sqrt.(x .^ 2 + y .^ 2)
    theta = atan.(y,x)
    dtheta_dy = cos.(theta) ./ r
    dr_dy = sin.(theta)
    dv_dr = 2 .* (r .- 1) .* cos.(theta) .+ sin.(theta)
    dv_dtheta = -1 .*  (r .- 1) .^ 2 .* sin.(theta) .+ (r .- 1) .* cos.(theta)
    return dv_dr .* dr_dy + dv_dtheta .* dtheta_dy
  end


  polar_laplace(x,y,e) = begin #u_rr + (1/r)*u_r + (1/r^2)*u_theta,theta
    r = sqrt.(x .^ 2 + y .^ 2)
    theta = atan.(y,x)
    if EToDomain[e] == 1
      u_r = c .* (2 .* r .^ 2 .* exp.(-1 .* r .^ 2) .+ 1 .- exp.(-1 .* r .^ 2)) .* sin.(theta)
      u_rr = c .* exp.(-1 .* r .^ 2) .* (6 .* r .- 4 .* r .^ 3) .* sin.(theta)
      return u_rr .+ (1 ./ r) .* u_r .- (c ./ r .^ 2) .* (1 .- exp.(-1 .* r .^ 2)) .* r .* sin.(theta)
      return dv_dr .* dr_dy .+ (cos.(r .- 1) .* cos.(theta)) .* (cos.(theta))
    elseif EToDomain[e] == 2
      return 2 * cos.(theta) .+ (1 ./ r) .* (2 .* (r .- 1) .* cos.(theta) .+ sin.(theta)) .+ (1 ./ r .^ 2) .* (-1 .* (r .- 1) .^ 2 .* cos.(theta) .- (r .- 1) .* sin.(theta))
    else
      error("invalid block")
    end
  end



  vex(x,y,e) = begin
    if EToDomain[e] == 1
      return vinside(x,y,e)
    elseif EToDomain[e] == 2
      return voutside(x,y,e)
    else
      error("invalid block")
    end
  end

  vex_x(x,y,e) = begin
    if EToDomain[e] == 1
      return vinside_x(x,y,e)
    elseif EToDomain[e] == 2
      return voutside_x(x,y,e)
    else
      error("invalid block")
    end
  end

  vex_y(x,y,e) = begin
    if EToDomain[e] == 1
      return vinside_y(x,y,e)
    elseif EToDomain[e] == 2
      return voutside_y(x,y,e)
    else
      error("invalid block")
    end
  end


  ϵ = zeros(4)

  for lvl = 1:length(ϵ)
    # Set up the local grid dimensions
    Nr = EToN0[1, :] * (2^(lvl-1))
    Ns = EToN0[2, :] * (2^(lvl-1))

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
    lvl == 1 && plot_blocks(lop)

    #
    # Do some assemble of the global volume operators
    #
    (M, FbarT, D, vstarts, FToλstarts) =
    LocalGlobalOperators(lop, Nr, Ns, FToB, FToE, FToLF, EToO, EToS,
                         (x) -> cholesky(Symmetric(x)))
    @show lvl
    locfactors = M.F

    # Get a unique array indexes for the face to jumps map
    FToδstarts = bcstarts(FToB, FToE, FToLF, BC_JUMP_INTERFACE, Nr, Ns)

    # Compute the number of volume, trace (λ), and jump (δ) points
    VNp = vstarts[nelems+1]-1
    λNp = FToλstarts[nfaces+1]-1
    δNp = FToδstarts[nfaces+1]-1
    # @show (VNp, λNp, δNp)

    # Build the (sparse) λ matrix using the schur complement and factor
    B = assembleλmatrix(FToλstarts, vstarts, EToF, FToB, locfactors, D, FbarT)
    BF = cholesky(Symmetric(B))

    (bλ, λ, gδ) = (zeros(λNp), zeros(λNp), zeros(λNp))
    (Δ, u, g) = (zeros(VNp), zeros(VNp), zeros(VNp))
    δ = zeros(δNp)


    for f = 1:nfaces
      if FToB[f] == BC_JUMP_INTERFACE
        (e1, e2) = FToE[:, f]
        (lf1, lf2) = FToLF[:, f]
        (xf, yf) = lop[e1].facecoord
        @views δ[FToδstarts[f]:(FToδstarts[f+1]-1)] =
        vex(xf[lf1], yf[lf1], e2) - vex(xf[lf1], yf[lf1], e1)
        #@show δ[FToδstarts[f]:(FToδstarts[f+1]-1)]
      end
    end

    bc_Dirichlet = (lf, x, y, e, δ) -> vex(x, y, e)
    bc_Neumann   = (lf, x, y, nx, ny, e, δ) -> (nx .* vex_x(x, y, e)
                                                + ny .* vex_y(x, y, e))
    in_jump      = (lf, x, y, e, δ) -> begin
      f = EToF[lf, e]
      if EToS[lf, e] == 1
        if EToO[lf, e]
          return -δ[FToδstarts[f]:(FToδstarts[f+1]-1)]
        else
          error("shouldn't get here")
        end
      else
        if EToO[lf, e]
          return  δ[FToδstarts[f]:(FToδstarts[f+1]-1)]
        else
          return  δ[(FToδstarts[f+1]-1):-1:FToδstarts[f]]
        end
      end
    end

    for e = 1:nelems
      gδe = ntuple(4) do lf
        f = EToF[lf, e]
        if EToO[lf, e]
          return @view gδ[FToλstarts[f]:(FToλstarts[f+1]-1)]
        else
          return  @view gδ[(FToλstarts[f+1]-1):-1:FToλstarts[f]]
        end
      end
      locbcarray!((@view g[vstarts[e]:vstarts[e+1]-1]), gδe, lop[e],
                  FToB[EToF[:,e]], bc_Dirichlet, bc_Neumann, in_jump, (e, δ))

      source = (x, y, e) -> (-polar_laplace(x, y, e))
      locsourcearray!((@view g[vstarts[e]:vstarts[e+1]-1]), source, lop[e], e)
    end

    if method == :monolithic
      # Solve the monolithic sytem
      M = blockdiag(ntuple(i->lop[i].M̃, length(lop))...)
      A = [M FbarT'; FbarT Diagonal(D)]
      uλ = A \ [g;gδ]
      u[:] .= uλ[1:length(u)]
    elseif method == :trace
      # Solve for the trace variables then compute displacements
      LocalToGLobalRHS!(bλ, g, gδ, u, locfactors, FbarT, vstarts)
      λ[:] = BF \ bλ

      u[:] = -FbarT' * λ
      u[:] .= g .+ u
      for e = 1:nelems
        F = locfactors[e]
        @views u[vstarts[e]:(vstarts[e+1]-1)] = F \ u[vstarts[e]:(vstarts[e+1]-1)]
        #=
        ldiv!((@view u[vstarts[e]:(vstarts[e+1]-1)]), F,
              (@view u[vstarts[e]:(vstarts[e+1]-1)]))
        =#
      end
    elseif method == :displacements
      # Solve for the discplacement directly
      M = blockdiag(ntuple(i->lop[i].M̃, length(lop))...)
      C = M - FbarT' * Diagonal(1 ./ D) * FbarT
      u[:] = C \ (g - FbarT' * (gδ ./ D))
    end

    for e = 1:nelems
      (x, y) = lop[e].coord
      JH = lop[e].JH
      #@show vex(x[:],y[:],e)' * JH * vex(x[:],y[:],e)
      @views Δ[vstarts[e]:(vstarts[e+1]-1)] = (u[vstarts[e]:(vstarts[e+1]-1)] -
                                               vex(x[:], y[:], e))
      ϵ[lvl] += (Δ[vstarts[e]:(vstarts[e+1]-1)]' * JH * Δ[vstarts[e]:(vstarts[e+1]-1)])

    end
    ϵ[lvl] = sqrt(ϵ[lvl])
    @show (lvl, ϵ[lvl])
  end
  println((log.(ϵ[1:end-1]) - log.(ϵ[2:end])) / log(2))

  nothing
end
