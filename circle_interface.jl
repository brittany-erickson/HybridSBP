include("global_curved.jl")

let
  # 3-----4-----4
  # |\         /|
  # | 11  2  12 |
  # |  \     /  |
  # |   7-8-8   |
  # |   |   |   |
  # 1 3 5 1 6 5 2
  # |   |   |   |
  # |   5-7-6   |
  # |  /     \  |
  # | 9   4  10 |
  # |/         \|
  # 1-----3-----2

  # Set up the connectivity
  # verts: Vertices
  verts = ((-1/√2, -1/√2), ( 1/√2, -1/√2), (-1/√2,  1/√2), ( 1/√2,  1/√2),
           (-1/3 , -1/3 ), ( 1/3 , -1/3 ), (-1/3 ,  1/3 ), ( 1/3 ,  1/3 ))

  # EToV: Element to Vertices
  EToV = ((5,6,7,8), (7,8,3,4), (1,5,3,7), (1,2,5,6), (6,2,8,4))

  # EToF: Element to Unique Global Faces
  EToF = ((5,6,7,8), (11,12,8,4), (1,5,9,11), (9,10,3,7), (6,2,10,12))

  # FToB: Unique Global Face to Boundary Conditions
  FToB = fill(BC_DIRICHLET, 12)
  FToB[5:8 ] .= BC_JUMP_INTERFACE
  FToB[9:12] .= BC_LOCKED_INTERFACE

  # EToN0: Element to base size sizes
  EToN0 = ((20, 20), (20, 12), (12, 20), (20, 12), (12, 20))

  # This is just needed because functions below expect arrays
  verts = flatten_tuples(verts)
  EToV  = flatten_tuples(EToV)
  EToF  = flatten_tuples(EToF)
  EToN0 = flatten_tuples(EToN0)

  # number of elements and faces
  (nelems, nfaces) = (size(EToV, 2), size(FToB, 1))

  # Some sanity checks
  @assert typeof(EToV) == Array{Int, 2} && size(EToV) == (4, nelems)
  @assert typeof(EToF) == Array{Int, 2} && size(EToF) == (4, nelems)
  @assert maximum(maximum(EToF)) == nfaces

  #{{{ Plot the connectivity using vertices corners
  @plotting (p1, p2, p3) = (plot(), plot(), plot())
  @plotting let
    # Do some plotting
    scatter!(p1, verts[1,:], verts[2,:], marker=10, legend=:none)
    for e = 1:nelems
      plot!(p1, verts[1, EToV[[1 2 4 3 1], e]]',
            verts[2, EToV[[1 2 4 3 1], e]]', legend=:none)
    end
  end
  #}}}

  # Determine secondary arrays
  # FToE : Unique Global Face to Element Number
  # FToLF: Unique Global Face to Element local face number
  # EToO : Element to Unique Global Faces Orientation
  # EToS : Element to Unique Global Face Side
  (FToE, FToLF, EToO, EToS) = connectivityarrays(EToV, EToF)

  # Exact solution
  (kx, ky) = (π, π)
  vex   = (x,y,e) ->       cos.(kx * x) .* cosh.(ky * y)
  vex_x = (x,y,e) -> -kx * sin.(kx * x) .* cosh.(ky * y)
  vex_y = (x,y,e) ->  ky * cos.(kx * x) .* sinh.(ky * y)
  if any(y -> y == BC_JUMP_INTERFACE, FToB)
    vex   = (x,y,e) ->       cos.(kx * x) .* cosh.(ky * y) .- 10*(e .== 1)
    vex_x = (x,y,e) -> -kx * sin.(kx * x) .* cosh.(ky * y)
    vex_y = (x,y,e) ->  ky * cos.(kx * x) .* sinh.(ky * y)
  end

  p = 4 # SBP interior order
  ϵ = zeros(2) # size of this array determines the number of levels to run

  OPTYPE = typeof(locoperator(2, 8, 8, (r,s)->r, (r,s)->s))
  for lvl = 1:length(ϵ)
    # Dictionary to store the operators
    lop = Dict{Int64, OPTYPE}()

    # Set up the local grid dimensions
    Nr = EToN0[1, :] * (2^(lvl-1))
    Ns = EToN0[2, :] * (2^(lvl-1))

    #{{{ Build the local volume operators
    for e = 1:nelems
      # Get the element corners
      (x1, x2, x3, x4) = verts[1, EToV[:, e]]
      (y1, y2, y3, y4) = verts[2, EToV[:, e]]

      if e == 1
        xt = (r,s)->transfinite_blend(x1, x2, x3, x4, r, s)
        yt = (r,s)->transfinite_blend(y1, y2, y3, y4, r, s)
      elseif e == 2
        xt = (r,s)->transfinite_blend(x1, x2, x3, x4, r, s;
                                      e4 = (α) -> cos.(π*(2 .- α)/4))
        yt = (r,s)->transfinite_blend(y1, y2, y3, y4, r, s;
                                      e4 = (α) -> sin.(π*(2 .- α)/4))
      elseif e == 3
        xt = (r,s)->transfinite_blend(x1, x2, x3, x4, r, s;
                                      e1 = (α) -> cos.(π*(4 .- α)/4))
        yt = (r,s)->transfinite_blend(y1, y2, y3, y4, r, s;
                                      e1 = (α) -> sin.(π*(4 .- α)/4))
      elseif e == 4
        xt = (r,s)->transfinite_blend(x1, x2, x3, x4, r, s;
                                      e3 = (α) -> cos.(π*(6 .+ α)/4))
        yt = (r,s)->transfinite_blend(y1, y2, y3, y4, r, s;
                                      e3 = (α) -> sin.(π*(6 .+ α)/4))
      elseif e == 5
        xt = (r,s)->transfinite_blend(x1, x2, x3, x4, r, s;
                                      e2 = (α) -> cos.(π*α/4))
        yt = (r,s)->transfinite_blend(y1, y2, y3, y4, r, s;
                                      e2 = (α) -> sin.(π*α/4))
      end

      # Build local operators
      lop[e] = locoperator(p, Nr[e], Ns[e], xt, yt, LFToB = FToB[EToF[:, e]])
    end
    #}}}

    # Assemble the global volume operators
    @plotting lvl == 1 && plotmesh(p2, lop, Nr, Ns, EToF, FToB)
    (vstarts, M, H, X, Y, E) = glovoloperator(lop, Nr, Ns)
    VNp = vstarts[nelems+1]-1

    # Build the trace operators
    (λstarts, T, D) = gloλoperator(lop, vstarts, FToB, FToE, FToLF, EToO, EToS,
                                   Nr, Ns)
    λNp = λstarts[nfaces+1]-1

    #{{{ Compute the boundary conditions
    bc_Dirichlet = (lf, x, y, e) -> vex(x, y, e)
    bc_Neumann   = (lf, x, y, nx, ny, e) -> (nx .* vex_x(x, y, e)
                                          + ny .* vex_y(x, y, e))
    in_jump      = (lf, x, y, e) -> begin
      f = EToF[lf, e]
      en = (EToS[lf, e] == 1 ? FToE[2,f] : FToE[1,f])
      @assert en != e
      (vex(x, y, e) - vex(x, y, en))
    end

    g = zeros(VNp)
    for e = 1:nelems
      locbcarray!((@view g[vstarts[e]:vstarts[e+1]-1]), lop[e], FToB[EToF[:,e]],
                  bc_Dirichlet, bc_Neumann, in_jump, (e,))
    end
    #}}}

    A = [ M -T'; -T  D ]
    uλ = A \ [g;zeros(λNp)]
    u = uλ[1:VNp]
    Δ = u - vex(X, Y, E)
    ϵ[lvl] = sqrt(Δ' * H * Δ)

    println("level = ", lvl, " :: error = ", ϵ[lvl])

    #{{{ Plot the solution
    @plotting if lvl == 1
      clims = (minimum(uλ[1:VNp]), maximum(uλ[1:VNp]))
      for e = 1:nelems
        (x, y) = lop[e][4]
        u = uλ[vstarts[e]:(vstarts[e+1]-1)]
        contour!(p3, reshape(x, Nr[e]+1, Ns[e]+1),
              reshape(y, Nr[e]+1, Ns[e]+1),
              reshape(u, Nr[e]+1, Ns[e]+1),
              fill=true, c = :balance, clims = clims,
              levels=linspace(clims[1], clims[2], 50))
      end
      plot!(p1, aspect_ratio = 1)
      plot!(p2, aspect_ratio = 1)
      plot!(p3, aspect_ratio = 1)
      display(plot(p1, p2, p3, layout = (1,3)))
    end
    #}}}
  end
  println((log.(ϵ[1:end-1]) - log.(ϵ[2:end])) / log(2))
end
