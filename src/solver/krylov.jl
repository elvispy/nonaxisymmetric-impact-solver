"""
    gmres_solve(system, b; x0=nothing)

Solve the matrix-free timestep system with GMRES. This is currently unpreconditioned.
"""
function gmres_solve(system::SolverSystem, b::Vector{Float64}; x0::Union{Nothing, Vector{Float64}}=nothing)
    cfg = system.config
    n = length(b)
    x = isnothing(x0) ? zeros(Float64, n) : copy(x0)

    r = b - apply_step_operator(system, x)
    β = norm(r)
    β < cfg.krylov_tol && return x

    maxiter = min(cfg.krylov_maxiter, n)
    V = zeros(Float64, n, maxiter + 1)
    H = zeros(Float64, maxiter + 1, maxiter)
    cs = zeros(Float64, maxiter)
    sn = zeros(Float64, maxiter)
    g = zeros(Float64, maxiter + 1)

    V[:, 1] .= r ./ β
    g[1] = β

    for j in 1:maxiter
        w = apply_step_operator(system, view(V, :, j))
        for i in 1:j
            H[i, j] = dot(view(V, :, i), w)
            w .-= H[i, j] .* view(V, :, i)
        end
        H[j + 1, j] = norm(w)
        if H[j + 1, j] != 0.0
            V[:, j + 1] .= w ./ H[j + 1, j]
        end

        for i in 1:(j - 1)
            temp = cs[i] * H[i, j] + sn[i] * H[i + 1, j]
            H[i + 1, j] = -sn[i] * H[i, j] + cs[i] * H[i + 1, j]
            H[i, j] = temp
        end

        ρ = hypot(H[j, j], H[j + 1, j])
        if ρ == 0.0
            cs[j] = 1.0
            sn[j] = 0.0
        else
            cs[j] = H[j, j] / ρ
            sn[j] = H[j + 1, j] / ρ
        end

        H[j, j] = cs[j] * H[j, j] + sn[j] * H[j + 1, j]
        H[j + 1, j] = 0.0

        temp = cs[j] * g[j] + sn[j] * g[j + 1]
        g[j + 1] = -sn[j] * g[j] + cs[j] * g[j + 1]
        g[j] = temp

        if abs(g[j + 1]) < cfg.krylov_tol
            y = H[1:j, 1:j] \ g[1:j]
            x .+= V[:, 1:j] * y
            return x
        end
    end

    ysmall = H[1:maxiter, 1:maxiter] \ g[1:maxiter]
    x .+= V[:, 1:maxiter] * ysmall
    return x
end
