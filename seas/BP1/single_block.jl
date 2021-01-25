include("../../global_curved.jl")
let
    # SBP interior order
    SBPp   = 2


    # mesh file side set type to actually boundary condition type
    bc_map = [BC_DIRICHLET, BC_DIRICHLET, BC_NEUMANN, BC_NEUMANN,
            BC_JUMP_INTERFACE]

    (verts, EToV, EToF, FToB, EToDomain) = read_inp_2d("meshes/1_1_block.inp";
            bc_map=bc_map)

    # number of elements and faces
    (nelems, nfaces) = (size(EToV, 2), size(FToB, 1))
    @show (nelems, nfaces)


    # This is the base mesh size in each dimension
    N1 = N0 = 8

     # EToN0 is the base mesh size (e.g., before refinement)
    EToN0 = zeros(Int64, 2, nelems)
    EToN0[1, :] .= N0
    EToN0[2, :] .= N1

    @assert typeof(EToV) == Array{Int,2} && size(EToV) == (4, nelems)
    @assert typeof(EToF) == Array{Int,2} && size(EToF) == (4, nelems)
    @assert maximum(maximum(EToF)) == nfaces


    # Exact solution
    Lx = 80
    Ly = 80

    vex(x, y, e) = begin
        if EToDomain[e] == 1
            return 1 .+ 0 .* x + 0 .* y
        end
    end

    vex_x(x, y, e) = begin
        if EToDomain[e] == 1
            return 0 .* x + 0 .* y
        end
    end

    vex_y(x, y, e) = begin
        if EToDomain[e] == 1
            return 0 .* x + 0 .* y
        end
    end


        # Set up the local grid dimensions
        Nr = EToN0[1, :] * (2^0)
        Ns = EToN0[2, :] * (2^0)


        # Dictionary to store the operators
        OPTYPE = typeof(locoperator(2, 8, 8))
        lop = Dict{Int64,OPTYPE}() # Be extra careful about the () here


        el_x = 10 # set it to be infinity to have even spread
        el_y = 10
        xt = (r,s) -> (el_x .* tan.(atan((Lx )/el_x).* (0.5*r .+ 0.5))  , el_x .* sec.(atan((Lx )/el_x).* (0.5*r .+ 0.5)).^2 * atan((Lx)/el_x) * 0.5 ,zeros(size(s)))
        yt = (r,s) -> (el_y .* tan.(atan((Ly )/el_y).* (0.5*s .+ 0.5))  , zeros(size(r)), el_y .* sec.(atan((Ly )/el_y).*(0.5*s .+ 0.5)) .^2 * atan((Ly )/el_y) * 0.5 )


        # create metrics
        metrics = create_metrics(SBPp, Nr[1], Ns[1], xt, yt) # not quite sure about this part
        # println("Metrics created successfully")

        # create local operator
        LFtoB = [BC_DIRICHLET,BC_DIRICHLET,BC_NEUMANN,BC_NEUMANN]
        lop = Dict{Int64,OPTYPE}() # This step to create a dict is essential
        lop[1] = locoperator(SBPp, Nr[1], Ns[1], metrics, LFtoB) # this function might not be correct

        # obtain M
        factorization = (x) -> cholesky(Symmetric(x))
        M = SBPLocalOperator1(lop, Nr[1], Ns[1], factorization)

        # obtain ge that stores boundary data
        ge = zeros((Nr[1]+1) * (Ns[1]+1))
        e = 1


        # boundary conditions
        bc_Dirichlet = (lf, x, y, e) -> vex(x, y, e)
        bc_Neumann   = (lf, x, y, nx, ny, e) -> (nx .* vex_x(x, y, e)
                                                + ny .* vex_y(x, y, e))

        # set right-hand side and solve
        locbcarray_mod!(ge, lop[e], LFtoB, bc_Dirichlet, bc_Neumann,(e))
        @show numerical_solution = M.F[e] \ ge
        

end
