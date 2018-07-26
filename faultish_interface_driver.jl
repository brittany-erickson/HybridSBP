do_plotting = true
include("global_curved.jl")
if VERSION <= v"0.6.999999"
  ldiv! = A_ldiv_B!
  cholesky = cholfact
  using CholmodSolve2
end

let

  # This is just needed because functions below expect arrays
  (verts, EToV, EToF, FToB) = read_inp_2d("meshes/BP1_V1.inp")

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
    scatter!(p1, verts[1,:], verts[2,:], marker=1, legend=:none)
    for e = 1:nelems
      plot!(p1, verts[1, EToV[[1 2 4 3 1], e]]',
            verts[2, EToV[[1 2 4 3 1], e]]', legend=:none)
    end
    plot!(p1, aspect_ratio = 1)
  end
  #}}}

  EToN0 = fill(13, (2, nelems))

  # Determine secondary arrays
  # FToE : Unique Global Face to Element Number
  # FToLF: Unique Global Face to Element local face number
  # EToO : Element to Unique Global Faces Orientation
  # EToS : Element to Unique Global Face Side
  (FToE, FToLF, EToO, EToS) = connectivityarrays(EToV, EToF)

  p   = 4 # SBP interior order
  lvl = 1 # Refinement

  # Dictionary to store the operators
  OPTYPE = typeof(locoperator(2, 8, 8, (r,s)->r, (r,s)->s))
  lop = Dict{Int64, OPTYPE}()

  # Set up the local grid dimensions
  Nr = EToN0[1, :] * (2^(lvl-1))
  Ns = EToN0[2, :] * (2^(lvl-1))

  #{{{ Build the local volume operators
  for e = 1:nelems
    # Get the element corners
    (x1, x2, x3, x4) = verts[1, EToV[:, e]]
    (y1, y2, y3, y4) = verts[2, EToV[:, e]]

    xt = (r,s)->transfinite_blend(x1, x2, x3, x4, r, s)
    yt = (r,s)->transfinite_blend(y1, y2, y3, y4, r, s)

    # Build local operators
    lop[e] = locoperator(p, Nr[e], Ns[e], xt, yt, LFToB = FToB[EToF[:, e]])
  end
  #}}}

  # Assemble the global volume operators
  @plotting lvl == 1 && let
    plotmesh(p2, lop, Nr, Ns, EToF, FToB)
    plot!(p2, aspect_ratio = 1)
  end

  (M, T, D, vstarts, λstarts) =
  LocalGlobalOperators(lop, Nr, Ns, FToB, FToE, FToLF, EToO, EToS,
                       (x) -> cholesky(Symmetric(x)))

  VNp = vstarts[nelems+1]-1
  λNp = λstarts[nfaces+1]-1
  println((VNp, λNp))
  Ttranspose = T'


  # factorization = (x) -> lufact(x)
  Afun = (Aλ, λ, u, vu) -> begin
    mul!(u, Ttranspose, λ)
    #=
    for e = 1:nelems
      F = M.F[e]
      vu[e] .= F \ vu[e]
    end
    =#
    for e = 1:nelems
      F = M.F[e]
      ldiv!(vu[e], F, vu[e])
    end
    mul!(Aλ, T, u)
    map!((x,y,z) -> x * y - z, Aλ, D, λ, Aλ)
  end
  bfun = (b, g, u) -> begin
    for e = 1:nelems
      F = M.F[e]
      @views u[vstarts[e]:(vstarts[e+1]-1)] = F \ g[vstarts[e]:(vstarts[e+1]-1)]
    end
    mul!(b, T, u)
    b
  end

  u = zeros(VNp)
  vu = Array{typeof(view(u, 1:10))}(undef, nelems)
  for e = 1:nelems
    vu[e] = view(u, vstarts[e]:(vstarts[e+1]-1))
  end
  bλ = zeros(λNp)
  λ = zeros(λNp)

  #{{{ Compute the boundary conditions
  W = 40
  bc_Dirichlet = (lf, x, y, e, t) -> zeros(size(x))
  bc_Neumann   = (lf, x, y, nx, ny, e, t) -> zeros(size(x))
  in_jump      = (lf, x, y, e, t) -> (EToS[lf, e] == 1 ? -t : t) * (y/W .+ 1)
  bc_Dirichlet = (lf, x, y, e, t) -> x
  in_jump      = (lf, x, y, e, t) -> zeros(size(x))

  g = zeros(VNp)
  (s1, s2, s3, s4) = (similar(λ), similar(λ), similar(λ), similar(λ))

  fun(x,y) = Afun(x,y,u,vu)
  tol = 1e-8
  iter = cg!(λ, s1, s2, s3, s4, bλ, (x,y)->Afun(x, y, u, vu); MaxIter=λNp,
             tol=tol)
  iter = cg!(λ, s1, s2, s3, s4, bλ, fun; MaxIter=λNp, tol=tol)
  λ[:] .= 0
  for t = linspace(0,1,10)
    for e = 1:nelems
      locbcarray!((@view g[vstarts[e]:vstarts[e+1]-1]), lop[e], FToB[EToF[:,e]],
                  bc_Dirichlet, bc_Neumann, in_jump, (e,t))
    end
    bfun(bλ, g, u)
    @time iter = cg!(λ, s1, s2, s3, s4, bλ, fun; MaxIter=λNp, tol = tol)
    println((iter, λNp))
    if iter < 0
      println("CG did not converge")
    end
  end

  u[:] = T' * λ
  u[:] .= g .+ u
  for e = 1:nelems
    F = M.F[e]
    @views u[vstarts[e]:(vstarts[e+1]-1)] = F \ u[vstarts[e]:(vstarts[e+1]-1)]
  end
  #}}}

  @plotting lvl == 1 && let
    clims = (minimum(u), maximum(u))
    for e = 1:nelems
      (x, y) = lop[e][4]
      up = u[vstarts[e]:(vstarts[e+1]-1)]
      plot!(p3, reshape(x, Nr[e]+1, Ns[e]+1),
            reshape(y, Nr[e]+1, Ns[e]+1),
            reshape(up, Nr[e]+1, Ns[e]+1),
            st = :surface, c = :balance, clims = clims)
    end
    # plot!(p3, aspect_ratio = 1, camera = (0, 90))
    # plot!(p3, aspect_ratio = 1)
    plot!(p3, aspect_ratio = 1, camera = (45, 45))

    @plotting display(plot(p1, p2, p3, layout = (3,1)))
  end
end
