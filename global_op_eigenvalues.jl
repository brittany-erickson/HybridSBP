# Force the random number generator to be the same for all runs of the script
using Random
Random.seed!(777)

# include("global_curved.jl")

import PGFPlots

let
  τs = 1

  verts = [-1 0 1 -1 0 1;
            0 0 0  1 1 1]
  EToV = [1 2 4 5;
          2 3 5 6]'
  EToF = [1 2 3 4;
          2 5 6 7]'
  FToB = fill(BC_DIRICHLET, (maximum(EToF),))
  FToB[2] = BC_LOCKED_INTERFACE

  (FToE, FToLF, EToO, EToS) = connectivityarrays(EToV, EToF)
  # plot_connectivity(verts, EToV)

  (nelems, nfaces) = (size(EToV, 2), size(FToB, 1))

  # Size of the grid
  SBP_orders = (2, 4, 6)

  num_samp = 1000
  min_eig_noSchur = zeros(length(SBP_orders), num_samp)
  max_eig_noSchur = zeros(length(SBP_orders), num_samp)
  min_eig_Schur_D = zeros(length(SBP_orders), num_samp)
  max_eig_Schur_D = zeros(length(SBP_orders), num_samp)
  min_eig_Schur_M = zeros(length(SBP_orders), num_samp)
  max_eig_Schur_M = zeros(length(SBP_orders), num_samp)
  for k = 1:num_samp
    k % 10 == 1 && println("sample $k of $num_samp")
    for (i, SBPp) in enumerate(SBP_orders)
      Ns = Nr = ntuple(e -> 3*SBPp - 1, nelems)

      # Build the local operators
      lop = ntuple(nelems) do e
        metrics = create_metrics(SBPp, Nr[e], Ns[e])

        # Build coefficients which are positive definite
        λ1 = rand(Nr[e]+1, Ns[e]+1)
        λ2 = rand(Nr[e]+1, Ns[e]+1) / 10000

        q  = π * rand(Nr[e]+1, Ns[e]+1)
        @. metrics.crr = λ1 * cos(q)^2 + λ2 * sin(q)^2
        @. metrics.css = λ1 * sin(q)^2 + λ2 * cos(q)^2
        @. metrics.crs = (λ2 - λ1) * cos(q) * sin(q)

        # Build local operators
        locoperator(SBPp, Nr[e], Ns[e], metrics; τscale = τs)
      end

      #
      # Assemble the global volume operators
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

      Np = (Nr.+1) .* (Ns.+1)
      A11 = [lop[1].M̃ spzeros(Np[1], Np[2])
             spzeros(Np[2], Np[1]) lop[2].M̃]
      A = [A11 FbarT';
           FbarT Diagonal(D)]


      # @show (VNp, λNp, δNp)
      # Build the (sparse) λ matrix using the schur complement and factor
      B = assembleλmatrix(FToλstarts, vstarts, EToF, FToB, locfactors, D, FbarT)
      @assert B ≈ Diagonal(D) - FbarT * (A11 \ Matrix(FbarT'))
      C = A11 - FbarT' * Diagonal(1 ./ D) * FbarT

      evA = eigen(Matrix(A)).values
      evB = eigen(Matrix(B)).values
      evC = eigen(Matrix(C)).values

      min_eig_noSchur[i, k], max_eig_noSchur[i, k] = extrema(evA)
      min_eig_Schur_M[i, k], max_eig_Schur_M[i, k] = extrema(evB)
      min_eig_Schur_D[i, k], max_eig_Schur_D[i, k] = extrema(evC)
    end
  end

  # Plot No Schur eigenvalues
  Lx = (1, num_samp)
  Ly = extrema(min_eig_noSchur)
  Ly = (floor(10000 * Ly[1]) / 10000, ceil(10000 * Ly[2]) / 10000)
  plt = Plot(BrailleCanvas(40, 20,
                           origin_x = Lx[1], origin_y = Ly[1], 
                           width = Lx[2] - Lx[1], height = Ly[2] - Ly[1]))

  annotate!(plt, :l, nrows(plt.graphics), string(Ly[1]), color = :light_black)
  annotate!(plt, :l, 1, string(Ly[2]), color = :light_black)
  annotate!(plt, :bl, string(Lx[1]), color = :light_black)
  annotate!(plt, :br, string(Lx[2]), color = :light_black)
  for (i, SBPp) = enumerate(SBP_orders)
    lineplot!(plt, min_eig_noSchur[i, :])
  end
  display(plt)

  pgf_axis = PGFPlots.Axis(style="width=5cm, height=5cm", xlabel="realization",
                           ylabel=PGFPlots.L"$\min \lambda$", xmin = 1,
                           xmax = num_samp)
  for (i, SBPp) = enumerate(SBP_orders)
    push!(pgf_axis, PGFPlots.Linear(1:num_samp, min_eig_noSchur[i, :],
                                    style="no marks, thick",
                                    legendentry = "order = $SBPp"))
  end
  PGFPlots.save("eigenvalues_noSchur.tikz", pgf_axis)

  # Plot Schur of D block eigenvalues
  Lx = (1, num_samp)
  Ly = extrema(min_eig_Schur_D)
  Ly = (floor(10000 * Ly[1]) / 10000, ceil(10000 * Ly[2]) / 10000)
  plt = Plot(BrailleCanvas(40, 20,
                           origin_x = Lx[1], origin_y = Ly[1], 
                           width = Lx[2] - Lx[1], height = Ly[2] - Ly[1]))

  annotate!(plt, :l, nrows(plt.graphics), string(Ly[1]), color = :light_black)
  annotate!(plt, :l, 1, string(Ly[2]), color = :light_black)
  annotate!(plt, :bl, string(Lx[1]), color = :light_black)
  annotate!(plt, :br, string(Lx[2]), color = :light_black)
  for (i, SBPp) = enumerate(SBP_orders)
    lineplot!(plt, min_eig_Schur_D[i, :])
  end
  display(plt)

  pgf_axis = PGFPlots.Axis(style="width=5cm, height=5cm", xlabel="realization",
                           ylabel=PGFPlots.L"$\min \lambda$", xmin = 1,
                           xmax = num_samp)
  for (i, SBPp) = enumerate(SBP_orders)
    push!(pgf_axis, PGFPlots.Linear(1:num_samp, min_eig_Schur_D[i, :],
                                    style="no marks, thick",
                                    legendentry = "order = $SBPp"))
  end
  PGFPlots.save("eigenvalues_Schur_D.tikz", pgf_axis)

  # Plot Schur of M block eigenvalues
  Lx = (1, num_samp)
  Ly = extrema(min_eig_Schur_M)
  Ly = (floor(10000 * Ly[1]) / 10000, ceil(10000 * Ly[2]) / 10000)
  plt = Plot(BrailleCanvas(40, 20,
                           origin_x = Lx[1], origin_y = Ly[1], 
                           width = Lx[2] - Lx[1], height = Ly[2] - Ly[1]))

  annotate!(plt, :l, nrows(plt.graphics), string(Ly[1]), color = :light_black)
  annotate!(plt, :l, 1, string(Ly[2]), color = :light_black)
  annotate!(plt, :bl, string(Lx[1]), color = :light_black)
  annotate!(plt, :br, string(Lx[2]), color = :light_black)
  for (i, SBPp) = enumerate(SBP_orders)
    lineplot!(plt, min_eig_Schur_M[i, :])
  end
  display(plt)

  pgf_axis = PGFPlots.Axis(style="width=5cm, height=5cm", xlabel="realization",
                           ylabel=PGFPlots.L"$\min \lambda$", xmin = 1,
                           xmax = num_samp)
  for (i, SBPp) = enumerate(SBP_orders)
    push!(pgf_axis, PGFPlots.Linear(1:num_samp, min_eig_Schur_M[i, :],
                                    style="no marks, thick",
                                    legendentry = "order = $SBPp"))
  end
  PGFPlots.save("eigenvalues_Schur_M.tikz", pgf_axis)
end
