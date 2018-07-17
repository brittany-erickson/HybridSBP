do_plotting = false
include("global_curved.jl")
if VERSION <= v"0.6.999999"
  ldiv! = A_ldiv_B!
  cholesky = cholfact
  # using CholmodSolve2
end

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

  # Set up the connectivity
  # verts: Vertices
  verts = ((0,1), (-1/2, 1/2), (1/2, 1/2),
           (0, 1/3), (-1, 0), (0,0), (1,0))

  # EToV: Element to Vertices
  EToV = ((1, 2, 3, 4), (5, 6, 2, 4), (6, 7, 4, 3))

  # EToF: Element to Unique Global Faces
  EToF = ((6, 7, 1, 9), (2, 8, 3, 7), (8, 5, 4, 9))

  # FToB: Unique Global Face to Boundary Conditions
  FToB = fill(BC_DIRICHLET, 9)
  FToB[2:7] .= BC_NEUMANN
  FToB[8] = BC_LOCKED_INTERFACE
  FToB[9] = BC_JUMP_INTERFACE

  # EToN0: Element to base size sizes
  EToN0 = ((16, 13), (14, 17), (16, 17))
  # EToN0 = ((15, 15), (15, 15), (15, 15))

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
    vex   = (x,y,e) ->       cos.(kx * x) .* cosh.(ky * y) .- 4*div.(e,2)
    vex_x = (x,y,e) -> -kx * sin.(kx * x) .* cosh.(ky * y)
    vex_y = (x,y,e) ->  ky * cos.(kx * x) .* sinh.(ky * y)
  end

  p = 4 # SBP interior order
  ϵ = zeros(5) # size of this array determines the number of levels to run

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

      # This will create the "curved" triangle by first creating the
      # straight-sided then using a global wrapping
      rt = (r,s)->transfinite_blend(x1, x2, x3, x4, r, s)
      st = (r,s)->transfinite_blend(y1, y2, y3, y4, r, s)

      xg = (r, s)-> r + sin.(π * s) .* cos.(π * r) / 8
      yg = (r, s)-> s - cos.(π * s) .* sin.(π * r) / 8

      xt = (r,s)->xg(rt(r,s), st(r,s))
      yt = (r,s)->yg(rt(r,s), st(r,s))

      #=
      # If you just want straight-sided elements use the functions below
      xt = (r,s)->transfinite_blend(x1, x2, x3, x4, r, s)
      yt = (r,s)->transfinite_blend(y1, y2, y3, y4, r, s)
      =#

      # Build local operators
      lop[e] = locoperator(p, Nr[e], Ns[e], xt, yt, LFToB = FToB[EToF[:, e]])
    end
    #}}}

    # Assemble the global volume operators
    @plotting lvl == 1 && plotmesh(p2, lop, Nr, Ns, EToF, FToB)
    (vstarts, ~, H, X, Y, E) = glovoloperator(lop, Nr, Ns)
    VNp = vstarts[nelems+1]-1

    # Build the trace operators
    (λstarts, T, D) = gloλoperator(lop, vstarts, FToB, FToE, FToLF, EToO, EToS,
                                   Nr, Ns)
    Ttranspose = T'
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

    # factorization = (x) -> lufact(x)
    factorization = (x) -> cholesky(Symmetric(x))
    FTYPE = typeof(factorization(sparse([1],[1],[1.0])))
    factors = Array{FTYPE, 1}(undef, nelems)
    for e = 1:nelems
      factors[e] = factorization(lop[e][1])
    end
    Afun = (Aλ, λ, u) -> begin
      mul!(u, Ttranspose, λ)
      for e = 1:nelems
        F = factors[e]
        @views u[vstarts[e]:(vstarts[e+1]-1)] = F \ u[vstarts[e]:(vstarts[e+1]-1)]
      end
      mul!(Aλ, T, u)
      map!((x,y,z) -> x * y - z, Aλ, D, λ, Aλ)
    end
    bfun = (b, g, u) -> begin
      for e = 1:nelems
        F = factors[e]
        @views u[vstarts[e]:(vstarts[e+1]-1)] = F \ g[vstarts[e]:(vstarts[e+1]-1)]
      end
      mul!(b, T, u)
      b
    end
    ufun = (u, λ) -> begin
    end

    u = zeros(VNp)
    bλ = zeros(λNp)
    (λ, iter) = cg(zeros(λNp), bfun(bλ, g, u), (x,y) -> Afun(x, y, u); MaxIter=λNp, tol = 1e-10)
    @time (λ, iter) = cg(zeros(λNp), bfun(bλ, g, u), (x,y) -> Afun(x, y, u); MaxIter=λNp, tol = 1e-10)
    if iter < 0
      println("CG did not converge")
    end
    u[:] = T' * λ
    u[:] .= g .+ u
    for e = 1:nelems
      F = factors[e]
      @views u[vstarts[e]:(vstarts[e+1]-1)] = F \ u[vstarts[e]:(vstarts[e+1]-1)]
    end
    Δ = u - vex(X, Y, E)
    ϵ[lvl] = sqrt(sum(H .* Δ.^2))
    println("level = ", lvl, " :: error = ", ϵ[lvl])

    #{{{ Plot the solution
    @plotting if lvl == 1
      clims = (minimum(u), maximum(u))
      for e = 1:nelems
        (x, y) = lop[e][4]
        up = u[vstarts[e]:(vstarts[e+1]-1)]
        plot!(p3, reshape(x, Nr[e]+1, Ns[e]+1),
              reshape(y, Nr[e]+1, Ns[e]+1),
              reshape(up, Nr[e]+1, Ns[e]+1),
              st = :surface, c = :balance, clims = clims)
      end
      plot!(p1, aspect_ratio = 1)
      plot!(p2, aspect_ratio = 1)
      plot!(p3, aspect_ratio = 1, camera = (0, 90))
      display(plot(p1, p2, p3, layout = (1,3)))
    end
    #}}}
  end
  println((log.(ϵ[1:end-1]) - log.(ϵ[2:end])) / log(2))
end
