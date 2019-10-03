include("global_curved.jl")

# This checks that the local solve is positive definite with random variable
# coefficients
let
  # SBP interior order
  SBPp   = 6

  Nr = 21
  Ns = 13
  metrics = create_metrics(SBPp, Nr, Ns)
  crr = metrics.crr
  css = metrics.css
  crs = metrics.crs

  for k = 1:1000
    # Build coefficients which are positive definite
    λ1 = rand(Nr+1, Ns+1)
    λ2 = rand(Nr+1, Ns+1)
    q  = rand(Nr+1, Ns+1)
    @. crr = λ1 * cos(q)^2 + λ2 * sin(q)^2
    @. css = λ1 * sin(q)^2 + λ2 * cos(q)^2
    @. crs = (λ2 - λ1) * cos(q) * sin(q)

    lop = locoperator(SBPp, Nr, Ns, metrics)

    @show extrema(eigen(Matrix(lop.M̃)).values)
  end
end
