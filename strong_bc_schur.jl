include("diagonal_sbp_D2.jl")
using SparseArrays

function make_operators_dirchlet(p, Nx, Ny)
  (D2x_1d, BSx_1d, HIx_1d, Hx_1d, rx_1d) = diagonal_sbp_D2(p, Nx)
  (D2y_1d, BSy_1d, HIy_1d, Hy_1d, ry_1d) = diagonal_sbp_D2(p, Ny)

  Ax_1d = BSx_1d - Hx_1d * D2x_1d
  Ay_1d = BSy_1d - Hy_1d * D2y_1d

  Ax = kron(Hy_1d, Ax_1d)
  Ay = kron(Ay_1d, Hx_1d)

  H  = kron(Hy_1d, Hx_1d)

  E0x_1d = sparse([     1], [     1], [1.], Nx + 1, Nx + 1)
  ENx_1d = sparse([Nx + 1], [Nx + 1], [1.], Nx + 1, Nx + 1)
  E0y_1d = sparse([     1], [     1], [1.], Ny + 1, Ny + 1)
  ENy_1d = sparse([Ny + 1], [Ny + 1], [1.], Ny + 1, Ny + 1)

  Inx = sparse(2:Nx,1:Nx-1, ones(Int64, Nx-1), Nx+1, Nx-1)
  Iny = sparse(2:Ny,1:Ny-1, ones(Int64, Ny-1), Ny+1, Ny-1)
  In = kron(Iny, Inx)

  Nfp = 2 * Nx + 2 * Ny
  Nxp = Nx + 1
  Nyp = Ny + 1
  B2C = [1; Nxp; 1+Nxp*Ny; Nxp*Nyp;
         ((Nxp+1):Nxp:Nxp*(Ny-1)+1);
         ((Nxp+1):Nxp:Nxp*(Ny-1)+1) .+ Nx;
         (2:Nx);
         (2:Nx) .+ Nxp*Ny]
  Bn = sparse(B2C, 1:Nfp, ones(Int64, Nfp), Nxp * Nyp, Nfp)

  x = kron(ones(Int64, Ny+1), rx_1d)
  y = kron(ry_1d, ones(Int64, Nx+1))

  (Ax, Ay, H, In, Bn, x, y)
end


let
  (Ax, Ay, H, In, Bn, x, y) = make_operators_dirchlet(2,3,4);
end
