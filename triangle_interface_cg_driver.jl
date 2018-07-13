do_plotting = false
include("global_curved.jl")
if VERSION <= v"0.6.999999"
  ldiv! = A_ldiv_B!
  cholesky = cholfact
end

if VERSION >= v"0.6.999999"
  # Transpose type bwdsub
  function bwdsub!(x, AT::Union{Transpose{T,SparseMatrixCSC{T,Int64}},
                                Adjoint{T,SparseMatrixCSC{T,Int64}}},
                   b) where T <: Real
    copy!(x, b)
    bwdsub!(AT, x)
  end
  function bwdsub!(AT::Union{Transpose{T,SparseMatrixCSC{T,Int64}},
                             Adjoint{T,SparseMatrixCSC{T,Int64}}},
                   b) where T <: Real
    nzval = AT.parent.nzval
    colval = AT.parent.rowval
    rowptr = AT.parent.colptr

    for i = length(b):-1:1
      jstart = rowptr[i]
      jend   = rowptr[i + 1] - 1

      # lopp through the row and subtract off pieces
      j = jend
      while j > jstart
        if colval[j] > i
          b[i] -= b[colval[j]]*nzval[j]
          j -= 1
        else
          break
        end
      end
      @assert colval[j] == i
      b[i] = b[i]/nzval[j]
    end
    b
  end
end

function bwdsub!(x, A::SparseMatrixCSC, b)
  copy!(x, b)
  bwdsub!(A, x)
end
function bwdsub!(A::SparseMatrixCSC, b)
  nzval = A.nzval
  rowval = A.rowval
  colptr = A.colptr

  for j = length(b):-1:1
    istart = colptr[j]
    iend   = colptr[j + 1] - 1
    while istart <= iend && rowval[iend] > j
      iend -= 1
    end

    @assert rowval[iend] == j

    b[j] = bj = b[j]/nzval[iend]

    # update remaining part
    for i = istart:iend-1
      b[rowval[i]] -= bj*nzval[i]
    end
  end
  b
end

function fwdsub!(x, A::SparseMatrixCSC, b)
  copy!(x, b)
  fwdsub!(A, x)
end
function fwdsub!(A::SparseMatrixCSC, b)
  nzval = A.nzval
  rowval = A.rowval
  colptr = A.colptr

  for j = 1:length(b)
    istart = colptr[j]
    iend   = colptr[j + 1] - 1
    while istart <= iend && rowval[istart] < j
      istart += 1
    end

    @assert rowval[istart] == j

    b[j] = bj = b[j]/nzval[istart]

    # update remaining part
    for i = istart+1:iend
      b[rowval[i]] -= bj*nzval[i]
    end
  end
  b
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
    (vstarts, M, H, X, Y, E) = glovoloperator(lop, Nr, Ns)
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

    #=
    utmp = zeros(VNp)
    CHOLTYPE = typeof(cholesky(sparse([1],[1],[1.0])))
    VIEWTYPE = typeof(view(utmp, 1:2))
    factors = Array{Tuple{CHOLTYPE, VIEWTYPE}, 1}(undef, nelems)
    for e = 1:nelems
      G = cholesky(Symmetric(lop[e][1]))
      V = view(utmp, vstarts[e]:(vstarts[e+1]-1))
      factors[e] = (G, V)
    end
    λtmp = zeros(λNp)
    dd = Vector(diag(D))
    Afun = (Aλ, λ) -> begin
      mul!(utmp, Ttranspose, λ)
      # utmp = M \ utmp
      for e = 1:nelems
        G = factors[e][1]
        v = factors[e][2]
        v[:] = G \ v
      end
      # Aλ[:] = D * λ - T * utmp
      mul!(λtmp, T, utmp)
      map!((x,y,z) -> x * y - z, Aλ, dd, λ, λtmp)
    end
    =#

    #=
    utmp = zeros(VNp)
    G = cholesky(sparse([1],[1],[1.0]))
    if VERSION <= v"0.6.999999"
      LTYPE = typeof(sparse(G[:L]))
      P = G[:p]
      LTTYPE = typeof(transpose(sparse(G[:L])))
    else
      LTYPE = typeof(sparse(G.L))
      P = G.p
      LTTYPE = typeof(transpose(sparse(G.L)))
    end
    L_factor = Array{LTYPE, 1}(undef, nelems)
    LT_factor = Array{LTTYPE, 1}(undef, nelems)
    VIEWTYPE = typeof(view(utmp, P))
    vtmp = Array{VIEWTYPE, 1}(undef, nelems)
    for e = 1:nelems
      G = cholesky(Symmetric(lop[e][1]))
      if VERSION <= v"0.6.999999"
        L = sparse(G[:L])
        P = G[:p]
        LT = transpose(sparse(G[:L]))
      else
        L = sparse(G.L)
        P = G.p
        LT = transpose(L)
      end
      L_factor[e] = L
      LT_factor[e] = LT
      vtmp[e] = view(utmp, P.+(vstarts[e]-1))
    end
    λtmp = zeros(λNp)
    dd = Vector(diag(D))
    Afun = (Aλ, λ) -> begin
      mul!(utmp, Ttranspose, λ)
      for e = 1:nelems
        bwdsub!(LT_factor[e], fwdsub!(L_factor[e], vtmp[e]))
      end
      mul!(λtmp, T, utmp)
      map!((x,y,z) -> x * y - z, Aλ, dd, λ, λtmp)
    end
    =#

    utmp = zeros(VNp)
    utmp2 = zeros(VNp)
    # factorization = (x) -> lufact(x)
    factorization = (x) -> cholesky(Symmetric(x))
    FTYPE = typeof(factorization(sparse([1],[1],[1.0])))
    VIEWTYPE = typeof(view(utmp, 1:2))
    factors = Array{Tuple{FTYPE, VIEWTYPE, VIEWTYPE}, 1}(undef, nelems)
    for e = 1:nelems
      F = factorization(lop[e][1])
      V = view(utmp, vstarts[e]:(vstarts[e+1]-1))
      V2 = view(utmp2, vstarts[e]:(vstarts[e+1]-1))
      factors[e] = (F, V, V2)
    end
    λtmp = zeros(λNp)
    dd = Vector(diag(D))
    Afun = (Aλ, λ) -> begin
      mul!(utmp, Ttranspose, λ)
      for e = 1:nelems
        F = factors[e][1]
        v = factors[e][2]
        v2 = factors[e][3]
        v2[:] = F \ v
        # ldiv!(v2, F, v)
      end
      # Aλ[:] = D * λ - T * utmp
      mul!(λtmp, T, utmp2)
      map!((x,y,z) -> x * y - z, Aλ, dd, λ, λtmp)
      # Aλ[:] .= 1
    end

    #=
    utmp = zeros(VNp)
    CHOLTYPE = typeof(cholesky(sparse([1],[1],[1.0])))
    VIEWTYPE = typeof(view(utmp, 1:2))
    factors = Array{Tuple{CHOLTYPE, VIEWTYPE}, 1}(undef, nelems)
    for e = 1:nelems
      G = cholesky(Symmetric(lop[e][1]))
      V = view(utmp, vstarts[e]:(vstarts[e+1]-1))
      factors[e] = (G, V)
    end
    λtmp = zeros(λNp)
    dd = Vector(diag(D))
    Afun = (Aλ, λ) -> begin
      mul!(utmp, Ttranspose, λ)
      # utmp = M \ utmp
      for e = 1:nelems
        G = factors[e][1]
        v = factors[e][2]
        v[:] = G \ v
      end
      # Aλ[:] = D * λ - T * utmp
      mul!(λtmp, T, utmp)
      map!((x,y,z) -> x * y - z, Aλ, dd, λ, λtmp)
    end
    =#

    (λ, iter) = cg(zeros(λNp), T * (M \ g), Afun; MaxIter=λNp, tol = 1e-10)
    if iter < 0
      println("CG did not converge")
    end
    u = M \ (g + T' * λ)
    Δ = u - vex(X, Y, E)
    ϵ[lvl] = sqrt(Δ' * H * Δ)
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
