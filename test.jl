do_plotting = false
if VERSION < v"0.6.999999"
  macro isdefined(s::Symbol)
    return isdefined(s)
  end
end
if !@isdefined mtime_global_curved
  mtime_global_curved = 0
end
if mtime_global_curved < mtime("global_curved.jl")
  println("including global_curved")
  include("global_curved.jl")
  mtime_global_curved = mtime("global_curved.jl")
end
using BenchmarkTools
let
  Nr = 50
  p = 4
  Nrp = Nr+1
  V = rand(Nrp)
  @benchmark variable_diagonal_sbp_D2($p, $Nr, $V)
end
