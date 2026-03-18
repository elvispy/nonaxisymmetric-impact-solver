using LinearAlgebra

struct SolverConfig
    nx::Int
    ny::Int
    dx::Float64
    dt::Float64
    steps::Int
    Fr::Float64
    We::Float64
    Re::Float64
    body_damping::Float64
    contact_mask::BitMatrix
    mass_density::Matrix{Float64}
    external_load::Vector{Float64}
    near_radius::Int
    quadrature_order::Int
    linear_solver::Symbol
    krylov_tol::Float64
    krylov_maxiter::Int
    include_surface_tension_boundary::Bool
end

struct SolverDomain
    nx::Int
    ny::Int
    N::Int
    dx::Float64
    x::Matrix{Float64}
    y::Matrix{Float64}
    laplacian::Matrix{Float64}
    dtn_operator::CorrectedKernelDtN2D
    contact_mask::BitMatrix
    contact_indices::Vector{Int}
    free_indices::Vector{Int}
    eta_free_map::Matrix{Float64}
    q_map::Matrix{Float64}
    pressure_map::Matrix{Float64}
    plane_basis::Matrix{Float64}
    force_map::Matrix{Float64}
    surface_tension_map::Matrix{Float64}
    mass_matrix::Matrix{Float64}
    gravity_load::Vector{Float64}
end

struct SolverSystem
    config::SolverConfig
    domain::SolverDomain
end

Base.@kwdef struct SolverState
    time::Float64
    eta::Matrix{Float64}
    phi::Matrix{Float64}
    pressure::Matrix{Float64}
    q::Vector{Float64}
    v::Vector{Float64}
end

"""
    SolverConfig(; ...)

Create the solver configuration for a prescribed-contact rigid-plane impact model.
The contact mask and mass density must match the bath grid dimensions.
"""
function SolverConfig(;
    nx::Integer,
    ny::Integer,
    dx::Real,
    dt::Real,
    steps::Integer,
    Fr::Real,
    We::Real,
    Re::Real,
    body_damping::Real=0.0,
    contact_mask::AbstractMatrix{Bool},
    mass_density::AbstractMatrix{<:Real},
    external_load::AbstractVector{<:Real}=zeros(3),
    near_radius::Integer=2,
    quadrature_order::Integer=8,
    linear_solver::Symbol=:gmres,
    krylov_tol::Real=1e-10,
    krylov_maxiter::Integer=200,
    include_surface_tension_boundary::Bool=true,
)
    nx_i = Int(nx)
    ny_i = Int(ny)
    size(contact_mask) == (nx_i, ny_i) || throw(DimensionMismatch("contact_mask must match grid size"))
    size(mass_density) == (nx_i, ny_i) || throw(DimensionMismatch("mass_density must match grid size"))
    length(external_load) == 3 || throw(DimensionMismatch("external_load must have length 3"))

    return SolverConfig(
        nx_i,
        ny_i,
        Float64(dx),
        Float64(dt),
        Int(steps),
        Float64(Fr),
        Float64(We),
        Float64(Re),
        Float64(body_damping),
        BitMatrix(contact_mask),
        Float64.(mass_density),
        Float64.(collect(external_load)),
        Int(near_radius),
        Int(quadrature_order),
        linear_solver,
        Float64(krylov_tol),
        Int(krylov_maxiter),
        include_surface_tension_boundary,
    )
end

"""
    zero_state(system)

Return the zero initial condition associated with a solver system.
"""
function zero_state(system::SolverSystem)
    nx = system.config.nx
    ny = system.config.ny
    return SolverState(
        0.0,
        zeros(nx, ny),
        zeros(nx, ny),
        zeros(nx, ny),
        zeros(3),
        zeros(3),
    )
end
