"""
    build_solver(config)

Build the full solver system, including geometry, operators, and body couplings.
"""
function build_solver(config::SolverConfig)
    domain = build_solver_domain(config)
    return SolverSystem(config, domain)
end

"""
    build_solver_domain(config)

Construct the spatial data needed by the solver: Cartesian Laplacian, DtN operator,
contact/free maps, generalized body mass matrix, and contact-boundary surface-tension map.
"""
function build_solver_domain(config::SolverConfig)
    x, y = grid_coordinates(config.nx, config.ny, config.dx)
    x_full = repeat(x, 1, config.ny)
    y_full = repeat(y, config.nx, 1)

    contact_indices = findall(vec(config.contact_mask))
    isempty(contact_indices) && throw(ArgumentError("contact_mask must contain at least one contact point"))
    free_indices = findall(.!vec(config.contact_mask))

    plane_basis = _plane_basis(x_full, y_full, contact_indices)
    rank(plane_basis) == 3 || throw(ArgumentError("contact_mask must contain at least three non-collinear points"))

    eta_free_map = zeros(Float64, config.nx * config.ny, length(free_indices))
    for (col, idx) in enumerate(free_indices)
        eta_free_map[idx, col] = 1.0
    end

    q_map = zeros(Float64, config.nx * config.ny, 3)
    for (row, idx) in enumerate(contact_indices)
        q_map[idx, :] .= plane_basis[row, :]
    end

    pressure_map = zeros(Float64, config.nx * config.ny, length(contact_indices))
    for (col, idx) in enumerate(contact_indices)
        pressure_map[idx, col] = 1.0
    end

    laplacian = _build_cartesian_laplacian(config.nx, config.ny, config.dx)
    dtn_operator = build_corrected_kernel_dtn(
        config.nx,
        config.ny,
        config.dx;
        near_radius=config.near_radius,
        quadrature_order=config.quadrature_order,
    )

    mass_matrix = _build_generalized_mass_matrix(x_full, y_full, config.mass_density, config.dx)
    force_map = (config.dx^2) .* transpose(plane_basis)
    surface_tension_map = _build_surface_tension_boundary_map(config.contact_mask, x_full, y_full, config.dx)
    gravity_load = _build_gravity_load(config.mass_density, x_full, y_full, config.dx, config.Fr)

    return SolverDomain(
        config.nx,
        config.ny,
        config.nx * config.ny,
        config.dx,
        x_full,
        y_full,
        laplacian,
        dtn_operator,
        config.contact_mask,
        contact_indices,
        free_indices,
        eta_free_map,
        q_map,
        pressure_map,
        plane_basis,
        force_map,
        surface_tension_map,
        mass_matrix,
        gravity_load,
    )
end

"""
    build_dense_dtn_matrix(nx, ny, dx; near_radius=2, quadrature_order=8)

Assemble a dense DtN matrix for small-grid reference solves and consistency tests.
"""
function build_dense_dtn_matrix(nx::Integer, ny::Integer, dx::Real; near_radius::Integer=2, quadrature_order::Integer=8)
    op = build_corrected_kernel_dtn(nx, ny, dx; near_radius=near_radius, quadrature_order=quadrature_order)
    n = Int(nx) * Int(ny)
    A = zeros(Float64, n, n)
    for j in 1:n
        basis = zeros(Float64, Int(nx), Int(ny))
        basis[j] = 1.0
        A[:, j] .= vec(apply(op, basis))
    end
    return A
end

"""Build the rigid-plane basis restricted to the contact-mask nodes."""
function _plane_basis(x::Matrix{Float64}, y::Matrix{Float64}, contact_indices::Vector{Int})
    basis = zeros(Float64, length(contact_indices), 3)
    for (row, idx) in enumerate(contact_indices)
        basis[row, 1] = 1.0
        basis[row, 2] = x[idx]
        basis[row, 3] = y[idx]
    end
    return basis
end

"""Assemble the generalized rigid-body mass matrix from a grid density field."""
function _build_generalized_mass_matrix(x::Matrix{Float64}, y::Matrix{Float64}, density::Matrix{Float64}, dx::Float64)
    basis = zeros(Float64, 3, length(density))
    basis[1, :] .= 1.0
    basis[2, :] .= vec(x)
    basis[3, :] .= vec(y)
    weights = vec(density) .* dx^2
    M = zeros(Float64, 3, 3)
    for i in 1:length(weights)
        M .+= weights[i] .* (basis[:, i] * transpose(basis[:, i]))
    end
    return M
end

"""Compute the generalized gravity load associated with the body mass distribution."""
function _build_gravity_load(density::Matrix{Float64}, x::Matrix{Float64}, y::Matrix{Float64}, dx::Float64, Fr::Float64)
    if !isfinite(Fr) || Fr == 0.0
        return zeros(3)
    end

    basis = (
        ones(length(density)),
        vec(x),
        vec(y),
    )
    weights = vec(density) .* dx^2 ./ Fr
    return -[
        sum(weights .* basis[1]),
        sum(weights .* basis[2]),
        sum(weights .* basis[3]),
    ]
end

"""Assemble the linear contact-boundary surface-tension map for the rigid-plane DOFs."""
function _build_surface_tension_boundary_map(contact_mask::BitMatrix, x::Matrix{Float64}, y::Matrix{Float64}, dx::Float64)
    nx, ny = size(contact_mask)
    A = zeros(Float64, 3, nx * ny)
    lin = LinearIndices((nx, ny))

    directions = (
        (-1, 0, -1.0, 0.0),
        (1, 0, 1.0, 0.0),
        (0, -1, 0.0, -1.0),
        (0, 1, 0.0, 1.0),
    )

    for j in 1:ny, i in 1:nx
        contact_mask[i, j] || continue
        inside_idx = lin[i, j]
        for (di, dj, nxn, nyn) in directions
            ii = i + di
            jj = j + dj
            outside = !(1 <= ii <= nx && 1 <= jj <= ny) || !contact_mask[ii, jj]
            outside || continue

            # Each exposed contact edge contributes a linear normal-slope term to the
            # generalized force/moment balance through its edge midpoint.
            midpoint_x = x[inside_idx] + 0.5 * dx * nxn
            midpoint_y = y[inside_idx] + 0.5 * dx * nyn
            basis = (1.0, midpoint_x, midpoint_y)

            A[:, inside_idx] .+= (-1.0) .* collect(basis)
            if 1 <= ii <= nx && 1 <= jj <= ny
                outside_idx = lin[ii, jj]
                A[:, outside_idx] .+= collect(basis)
            end
        end
    end

    return A
end

"""Build the dense 5-point Cartesian Laplacian on the full bath grid."""
function _build_cartesian_laplacian(nx::Int, ny::Int, dx::Float64)
    n = nx * ny
    L = zeros(Float64, n, n)
    invdx2 = 1.0 / dx^2

    for j in 1:ny, i in 1:nx
        idx = LinearIndices((nx, ny))[i, j]
        L[idx, idx] = -4.0 * invdx2
        if i > 1
            L[idx, LinearIndices((nx, ny))[i - 1, j]] = invdx2
        end
        if i < nx
            L[idx, LinearIndices((nx, ny))[i + 1, j]] = invdx2
        end
        if j > 1
            L[idx, LinearIndices((nx, ny))[i, j - 1]] = invdx2
        end
        if j < ny
            L[idx, LinearIndices((nx, ny))[i, j + 1]] = invdx2
        end
    end

    return L
end
