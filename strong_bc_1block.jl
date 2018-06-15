include("diagonal_sbp_D2.jl")

function make_operators_dirchlet(p, Nx, Ny)
  (D2x_1d, BSx_1d, HIx_1d, Hx_1d, rx_1d) = diagonal_sbp_D2(p, Nx)
  (D2y_1d, BSy_1d, HIy_1d, Hy_1d, ry_1d) = diagonal_sbp_D2(p, Ny)

  Ax_1d = BSx_1d - Hx_1d * D2x_1d
  Ay_1d = BSy_1d - Hy_1d * D2y_1d

  Ax = kron(Hy_1d, Ax_1d)
  Ay = kron(Ay_1d, Hx_1d)

  E0x_1d = sparse([     1], [     1], [1.], Nx + 1, Nx + 1)
  ENx_1d = sparse([Nx + 1], [Nx + 1], [1.], Nx + 1, Nx + 1)
  E0y_1d = sparse([     1], [     1], [1.], Ny + 1, Ny + 1)
  ENy_1d = sparse([Ny + 1], [Ny + 1], [1.], Ny + 1, Ny + 1)

  Inx = sparse(2:Nx,1:Nx-1, ones(Nx-1), Nx+1, Nx-1)
  Iny = sparse(2:Ny,1:Ny-1, ones(Ny-1), Ny+1, Ny-1)
  In = kron(Iny, Inx)

  Bn = I - In * In'

  x = kron(ones(Ny+1), rx_1d)
  y = kron(ry_1d, ones(Nx+1))

  (Ax, Ay, In, Bn, x, y)
end

let
  Nx = 30
  Ny = 60
  (Ax, Ay, In, Bn, x, y) = make_operators_dirchlet(10, Nx, Ny)
  display(reshape(Bn * x, Nx+1, Ny+1))
  display(reshape(Bn * y, Nx+1, Ny+1))
  display(reshape(In' * x, Nx-1, Ny-1))
  display(reshape(In' * y, Nx-1, Ny-1))
end
