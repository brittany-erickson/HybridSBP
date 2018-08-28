using Plots
let
  p1 = plot()
  year_seconds = 31556926.
  #=
  for (base_name, title) = (("BP1_small_SBPp2_ptsc100", "many blocks"),
                            ("BP1_small_SBPp2_ptsc100_lvl2", "many blocks, refined"),
                            ("BP1_uniform_small4_SBPp2_ptsc100", "4 blocks"),
                            ("BP1_uniform_small4_SBPp2_ptsc200", "4 blocks, bigger penalty"),
                            ("BP1_uniform_small4_SBPp2_ptsc50", "4 blocks, smaller penalty"),
                            ("BP1_uniform_small_SBPp2_ptsc100", "2 blocks"),
                            ("BP1_uniform_small_SBPp2", "2 blocks, smaller penalty"),
                            ("BP1_uniform_small", "2 blocks, 4th order (?)"))
  =#
  #=
  for (base_name, title) = (("BP1_V0_nodeepload_p_2_lvl_2", "no deep load"),
                            ("BP1_V0_nodeepload_skew_p_2_lvl_1", "no deep load, skew"))
  =#
  for (base_name, title) = (
                            # ("BP1_V0_nodeepload_p_2_lvl_2", "no deep load"),
                            # ("BP1_uniform_small4_SBPp2_ptsc100_lvl2", "4 block, refined"),
                            ("BP1_uniform_small4_SBPp2_ptsc100_lvl1_Lx36", "Lx36"),
                            # ("BP1_uniform_small_SBPp2_ptsc100_lvl1_Lx36", "Lx36"),
                           )
    @show (base_name, title)
    y = open("$(base_name)_slip.dat") do f
      node_data = split(readline(f))
      Np = length(node_data) - 2
      y = zeros(Np)
      for d = 1:Np
        y[d] = parse(Float64, node_data[d+2])
      end
      y
    end
    (tδ, header) = readdlm("$(base_name)_slip.dat"; header=true)

    p1 = plot()
    tlast = 0
    for t = 2:size(tδ, 1)
      if tδ[t, 1] - tδ[t-1, 1] > 0.1 * year_seconds
        if tδ[t,1] - tlast > year_seconds
          plot!(p1, tδ[t, 3:end], y, color=:blue, legend=:none)
          tlast = tδ[t,1]
        end
      else
        if tδ[t,1] - tlast > 1
          plot!(p1, tδ[t, 3:end], y, color=:red, legend=:none)
          tlast = tδ[t,1]
        end
      end
    end
    display(plot!(p1, ylims=(-20, 0), xlims=(0, 20), title=title))
    savefig("$(base_name).pdf")
  end
end
