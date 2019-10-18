# Force the random number generator to be the same for all runs of the script
using Random
Random.seed!(777)

include("global_curved.jl")

import PGFPlots

# This checks that the local solve is positive definite with random variable
# coefficients
let
  τs = 1

  # Size of the grid
  Nr = Ns = 3 .* (3, 5, 7) .- 1
  SBP_orders = (2, 4, 6)

  # Create the random variable coefficients
  metrics = ntuple(i -> create_metrics(SBP_orders[i], Nr[i], Ns[i]),
                   length(SBP_orders))
  crr = ntuple(i->metrics[i].crr, length(SBP_orders))
  css = ntuple(i->metrics[i].css, length(SBP_orders))
  crs = ntuple(i->metrics[i].crs, length(SBP_orders))

  num_samp = 1000
  min_eig = zeros(length(SBP_orders), 2, num_samp)
  for k = 1:num_samp
    k % 10 == 1 && println("sample $k of $num_samp")

    for (i, SBPp) in enumerate(SBP_orders)
      # Build coefficients which are positive definite
      λ1 = rand(Nr[i]+1, Ns[i]+1)
      λ2 = rand(Nr[i]+1, Ns[i]+1) / 10000

      q  = π * rand(Nr[i]+1, Ns[i]+1)
      @. crr[i] = λ1 * cos(q)^2 + λ2 * sin(q)^2
      @. css[i] = λ1 * sin(q)^2 + λ2 * cos(q)^2
      @. crs[i] = (λ2 - λ1) * cos(q) * sin(q)

      lop = locoperator(SBPp, Nr[i], Ns[i], metrics[i]; τscale = τs)
      eigs = eigen(Matrix(lop.M̃))
      min_eig[i, 1, k] = minimum(eigs.values)
      min_eig[i, 1, k] < 0 && @show (1, k, SBPp, min_eig[i, 1, k])

      lop = locoperator(SBPp, Nr[i], Ns[i], metrics[i],
                        (BC_DIRICHLET, BC_NEUMANN, BC_NEUMANN, BC_NEUMANN);
                        τscale = τs)
      eigs = eigen(Matrix(lop.M̃))
      min_eig[i, 2, k] = minimum(eigs.values)
      min_eig[i, 2, k] < 0 && @show (2, k, SBPp, min_eig[i, 2, k])
    end
  end

  # Dirichlet
  Lx = (1, num_samp)
  Ly = extrema(min_eig[:, 1, :])
  Ly = (floor(10000 * Ly[1]) / 10000, ceil(10000 * Ly[2]) / 10000)
  plt = Plot(BrailleCanvas(40, 20,
                           origin_x = Lx[1], origin_y = Ly[1], 
                           width = Lx[2] - Lx[1], height = Ly[2] - Ly[1]))

  annotate!(plt, :l, nrows(plt.graphics), string(Ly[1]), color = :light_black)
  annotate!(plt, :l, 1, string(Ly[2]), color = :light_black)
  annotate!(plt, :bl, string(Lx[1]), color = :light_black)
  annotate!(plt, :br, string(Lx[2]), color = :light_black)
  for (i, SBPp) = enumerate(SBP_orders)
    lineplot!(plt, min_eig[i, 1, :])
  end
  display(plt)

  pgf_axis = PGFPlots.Axis(style="width=5cm, height=5cm", xlabel="realization",
                           ylabel=PGFPlots.L"$\min \lambda$", xmin = 1,
                           xmax = num_samp)
  for (i, SBPp) = enumerate(SBP_orders)
    push!(pgf_axis, PGFPlots.Linear(1:num_samp, min_eig[i, 1, :],
                                    style="no marks, thick",
                                    legendentry = "order = $SBPp"))
  end
  PGFPlots.save("local_eigenvalues_Dirichlet.tikz", pgf_axis)

  # Dirichlet & Neumann
  Lx = (1, num_samp)
  Ly = extrema(min_eig[:, 2, :])
  Ly = (floor(10000 * Ly[1]) / 10000, ceil(10000 * Ly[2]) / 10000)
  plt = Plot(BrailleCanvas(40, 20,
                           origin_x = Lx[1], origin_y = Ly[1], 
                           width = Lx[2] - Lx[1], height = Ly[2] - Ly[1]))

  annotate!(plt, :l, nrows(plt.graphics), string(Ly[1]), color = :light_black)
  annotate!(plt, :l, 1, string(Ly[2]), color = :light_black)
  annotate!(plt, :bl, string(Lx[1]), color = :light_black)
  annotate!(plt, :br, string(Lx[2]), color = :light_black)
  for (i, SBPp) = enumerate(SBP_orders)
    lineplot!(plt, min_eig[i, 2, :])
  end
  display(plt)

  pgf_axis = PGFPlots.Axis(style="width=5cm, height=5cm", xlabel="realization",
                           ylabel=PGFPlots.L"$\min \lambda$", xmin = 1,
                           xmax = num_samp)
  for (i, SBPp) = enumerate(SBP_orders)
    push!(pgf_axis, PGFPlots.Linear(1:num_samp, min_eig[i, 2, :],
                                    style="no marks, thick",
                                    legendentry = "order = $SBPp"))
  end
  PGFPlots.save("local_eigenvalues_Neumann.tikz", pgf_axis)
end

# This checks that the local solve is positive definite with random variable
# coefficients
let
  # Size of the grid
  Nr = Ns = 3 .* (3, 5, 7) .- 1
  SBP_orders = (2, 4, 6)

  # Create the random variable coefficients
  metrics = ntuple(i -> create_metrics(SBP_orders[i], Nr[i], Ns[i]),
                   length(SBP_orders))
  crr = ntuple(i->metrics[i].crr, length(SBP_orders))
  css = ntuple(i->metrics[i].css, length(SBP_orders))
  crs = ntuple(i->metrics[i].crs, length(SBP_orders))
  metrics = ntuple(i -> create_metrics(SBP_orders[i], Nr[i], Ns[i]),
                   length(SBP_orders))
  crr = ntuple(i->metrics[i].crr, length(SBP_orders))
  css = ntuple(i->metrics[i].css, length(SBP_orders))
  crs = ntuple(i->metrics[i].crs, length(SBP_orders))

  for (i, SBPp) in enumerate(SBP_orders)
    # Build coefficients which are positive definite
    λ1 = rand(Nr[i]+1, Ns[i]+1)
    λ2 = rand(Nr[i]+1, Ns[i]+1)

    q  = rand(Nr[i]+1, Ns[i]+1)
    @. crr[i] = λ1 * cos(q)^2 + λ2 * sin(q)^2
    @. css[i] = λ1 * sin(q)^2 + λ2 * cos(q)^2
    @. crs[i] = (λ2 - λ1) * cos(q) * sin(q)
  end

  τscales = range(0.;stop=2, length=10)
  max_eig = zeros(length(SBP_orders), length(τscales))
  min_eig = zeros(length(SBP_orders), length(τscales))

  for (k, τs) in enumerate(τscales)
    k % 10 == 1 && println("sample $k of $(length(τscales))")

    for (i, SBPp) in enumerate(SBP_orders)
      lop = locoperator(SBPp, Nr[i], Ns[i], metrics[i]; τscale = τs)
      eigs = eigen(Matrix(lop.M̃))
      max_eig[i,k] = maximum(eigs.values)
      min_eig[i,k] = minimum(eigs.values)
    end
  end

  Lx = (τscales[1], τscales[end])
  Ly = extrema(real.(max_eig))
  Ly = (floor(Ly[1]), ceil(Ly[2]))
  plt_max = Plot(BrailleCanvas(40, 20,
                               origin_x = Lx[1], origin_y = Ly[1], 
                               width = Lx[2] - Lx[1], height = Ly[2] - Ly[1]))

  annotate!(plt_max, :l, nrows(plt_max.graphics), string(Ly[1]), color = :light_black)
  annotate!(plt_max, :l, 1, string(Ly[2]), color = :light_black)
  annotate!(plt_max, :bl, string(Lx[1]), color = :light_black)
  annotate!(plt_max, :br, string(Lx[2]), color = :light_black)
  for (i, SBPp) = enumerate(SBP_orders)
    lineplot!(plt_max, τscales, max_eig[i,:])
  end
  display(plt_max)

  pgf_axis = PGFPlots.Axis(style="width=5cm, height=5cm",
                           legendPos="north west",
                           xlabel=PGFPlots.L"$\tau_{s}$",
                           ylabel=PGFPlots.L"$\max \lambda$",
                           xmin = τscales[1],
                           xmax = τscales[end])
  for (i, SBPp) = enumerate(SBP_orders)
    push!(pgf_axis, PGFPlots.Linear(τscales, max_eig[i, :],
                                    style="no marks, thick",
                                    legendentry = "order = $SBPp"))
  end
  PGFPlots.save("tau_scaling_max_eig.tikz", pgf_axis)

  Lx = (τscales[1], τscales[end])
  Ly = extrema(min_eig)
  Ly = (floor(1000*Ly[1])/1000, ceil(1000*Ly[2])/1000)
  plt_max = Plot(BrailleCanvas(40, 20,
                               origin_x = Lx[1], origin_y = Ly[1], 
                               width = Lx[2] - Lx[1], height = Ly[2] - Ly[1]))

  annotate!(plt_max, :l, nrows(plt_max.graphics), string(Ly[1]), color = :light_black)
  annotate!(plt_max, :l, 1, string(Ly[2]), color = :light_black)
  annotate!(plt_max, :bl, string(Lx[1]), color = :light_black)
  annotate!(plt_max, :br, string(Lx[2]), color = :light_black)
  for (i, SBPp) = enumerate(SBP_orders)
    lineplot!(plt_max, τscales, min_eig[i,:])
  end
  display(plt_max)

  pgf_axis = PGFPlots.Axis(style="width=5cm, height=5cm",
                           legendPos="north west",
                           xlabel=PGFPlots.L"$\tau_{s}$",
                           ylabel=PGFPlots.L"$\max \lambda$",
                           xmin = τscales[1],
                           xmax = τscales[end])
  for (i, SBPp) = enumerate(SBP_orders)
    push!(pgf_axis, PGFPlots.Linear(τscales, min_eig[i, :],
                                    style="no marks, thick",
                                    legendentry = "order = $SBPp"))
  end
  PGFPlots.save("tau_scaling_min_eig.tikz", pgf_axis)
end
