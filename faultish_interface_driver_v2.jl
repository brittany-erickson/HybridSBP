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


  EToN0 = fill(13, (2, nelems))
  println(nelems)

  # Determine secondary arrays
  # FToE : Unique Global Face to Element Number
  # FToLF: Unique Global Face to Element local face number
  # EToO : Element to Unique Global Faces Orientation
  # EToS : Element to Unique Global Face Side
  (FToE, FToLF, EToO, EToS) = connectivityarrays(EToV, EToF)

  #{{{ Plot the connectivity using vertices corners
  @plotting (p1, p2, p3) = (plot(), plot(), plot())
  @plotting let
    # Do some plotting
    scatter!(p1, verts[1,:], verts[2,:], marker=1, legend=:none)
    LFToLV = flatten_tuples(((1,3), (2, 4), (1,2), (3,4)))
    for f = 1:nfaces
      e  = FToE[1,f]
      lf = FToLF[1,f]
      V = EToV[LFToLV[:,lf], e]
      if FToB[f] == BC_DIRICHLET
        plot!(p1, verts[1, V], verts[2, V], color=:red, linewidth=3)
      elseif FToB[f] == BC_NEUMANN
        plot!(p1, verts[1, V], verts[2, V], color=:blue, linewidth=3)
      elseif FToB[f] == BC_LOCKED_INTERFACE
        plot!(p1, verts[1, V], verts[2, V], color=:black, linewidth=1)
      elseif FToB[f] == BC_JUMP_INTERFACE
        plot!(p1, verts[1, V], verts[2, V], color=:green, linewidth=3)
      else
        error("invalid bc")
      end
    end
    plot!(p1, aspect_ratio = 1)
  end
  #}}}

  p   = 4 # SBP interior order
  lvl = 2 # Refinement

  # Dictionary to store the operators
  OPTYPE = typeof(locoperator(2, 8, 8, (r,s)->r, (r,s)->s))
  lop = Dict{Int64, OPTYPE}()

  # Set up the local grid dimensions
  Nr = EToN0[1, :] * (2^(lvl-1))
  Ns = EToN0[2, :] * (2^(lvl-1))

  #{{{ Build the local volume operators
  for e = 1:nelems
    # println((e, nelems))
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

  (M, T, D, vstarts, λstarts) =
  LocalGlobalOperators(lop, Nr, Ns, FToB, FToE, FToLF, EToO, EToS,
                       (x) -> cholesky(Symmetric(x)))

  jumpstarts = similar(λstarts)
  jumpstarts[1] = 1
  for f = 1:nfaces
    if FToB[f] == BC_JUMP_INTERFACE
      jumpstarts[f+1] = jumpstarts[f] + (λstarts[f+1]-λstarts[f])
    else
      jumpstarts[f+1] = jumpstarts[f]
    end
  end

  VNp = vstarts[nelems+1]-1
  λNp = λstarts[nfaces+1]-1
  jmpNp = jumpstarts[nfaces+1]-1
  println((VNp, λNp, jmpNp))
  Ttranspose = T'

  sz = λNp
  for e = 1:nelems
    lλs = Array{Int64, 1}(undef, 4)
    for lf = 1:4
      f = EToF[lf,e]
      lλs[lf] = λstarts[f+1] - λstarts[f]
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
  for e = 1:nelems
    # println((e, nelems))
    F = M.F[e]
    vrng = vstarts[e]:(vstarts[e+1]-1)
    @time for lf = 1:4
      f = EToF[lf,e]
      if FToB[f] == BC_LOCKED_INTERFACE || FToB[f] == BC_JUMP_INTERFACE
        λrng = λstarts[f]:(λstarts[f+1]-1)
        B = Matrix(F \ Ttranspose[vrng, λrng])
        for lf2 = 1:4
          f2 = EToF[lf2,e]
          if FToB[f2] == BC_LOCKED_INTERFACE || FToB[f2] == BC_JUMP_INTERFACE
            λrng2 = λstarts[f2]:(λstarts[f2+1]-1)
            C = T[λrng2, vrng] * B
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
  @time B = sparse(Ie, Je, Ve, λNp, λNp)
  @assert B ≈ B'
  println((λNp * λNp, nnz(B), nnz(B) / λNp^2))
  @time BF = cholesky(Symmetric(B))

  bfun = (b, g, u) -> begin
    for e = 1:nelems
      F = M.F[e]
      @views u[vstarts[e]:(vstarts[e+1]-1)] = F \ g[vstarts[e]:(vstarts[e+1]-1)]
    end
    mul!(b, T, u)
    b
  end

  bλ = zeros(λNp)
  λ = zeros(λNp)
  u = zeros(VNp)
  g = zeros(VNp)

  #{{{ Compute the boundary conditions
  δ = ones(jmpNp)
  W = 40
  bc_Dirichlet = (lf, x, y, e, t) -> zeros(size(x))
  bc_Neumann   = (lf, x, y, nx, ny, e, t) -> zeros(size(x))
  in_jump      = (lf, x, y, e, t) -> begin
    f = EToF[lf, e]
    if EToS[lf, e] == 1
      return δ[jumpstarts[f]:(jumpstarts[f+1]-1)]
    else
      return -δ[jumpstarts[f]:(jumpstarts[f+1]-1)]
    end
  end
  #}}}

  λ[:] .= 0
  # for t = linspace(0,1,10)
  for t = 1
    for e = 1:nelems
      locbcarray!((@view g[vstarts[e]:vstarts[e+1]-1]), lop[e], FToB[EToF[:,e]],
                  bc_Dirichlet, bc_Neumann, in_jump, (e,t))
    end
    λ[:] = BF \ bfun(bλ, g, u)
  end

  u[:] = T' * λ
  u[:] .= g .+ u
  for e = 1:nelems
    F = M.F[e]
    @views u[vstarts[e]:(vstarts[e+1]-1)] = F \ u[vstarts[e]:(vstarts[e+1]-1)]
  end

  @plotting let
    clims = (minimum(u), maximum(u))
    for e = 1:nelems
      (x, y) = lop[e][4]
      up = u[vstarts[e]:(vstarts[e+1]-1)]
      plot!(p2, reshape(x, Nr[e]+1, Ns[e]+1),
            reshape(y, Nr[e]+1, Ns[e]+1),
            reshape(up, Nr[e]+1, Ns[e]+1),
            st = :surface, c = :balance, clims = clims)
    end
    # plot!(p2, aspect_ratio = 1, camera = (0, 90))
    # plot!(p2, aspect_ratio = 1)
    plot!(p2, aspect_ratio = 1, camera = (45, 45))

    @plotting display(plot(p1, p2, layout = (2,1)))
  end
  nothing
end
