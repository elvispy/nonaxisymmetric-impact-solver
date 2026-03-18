using FFTW

const DTN_CONSTANT = 1.0 / (2.0 * π)

struct PeriodicDtN2D
    nx::Int
    ny::Int
    dx::Float64
    kernel::Matrix{Float64}
    kernel_fft::Matrix{ComplexF64}
end

struct CorrectedKernelDtN2D
    nx::Int
    ny::Int
    dx::Float64
    near_radius::Int
    quadrature_order::Int
    kernel::Matrix{Float64}
    kernel_fft::Matrix{ComplexF64}
    fft_shape::Tuple{Int, Int}
end

"""
    build_toeplitz_dtn(nx, ny, dx; near_radius=2, quadrature_order=8)

Build the default nonperiodic corrected-kernel DtN operator on a uniform 2D grid.
"""
function build_toeplitz_dtn(nx::Integer, ny::Integer, dx::Real; near_radius::Integer=2, quadrature_order::Integer=8)
    return build_corrected_kernel_dtn(nx, ny, dx; near_radius=near_radius, quadrature_order=quadrature_order)
end

"""
    build_corrected_kernel_dtn(nx, ny, dx; near_radius=2, quadrature_order=8)

Construct the corrected nonperiodic DtN operator, including the local singular correction
and the FFT embedding used for linear convolution.
"""
function build_corrected_kernel_dtn(nx::Integer, ny::Integer, dx::Real; near_radius::Integer=2, quadrature_order::Integer=8)
    nx_i = Int(nx)
    ny_i = Int(ny)
    dx_f = Float64(dx)
    near_i = Int(near_radius)
    quad_i = Int(quadrature_order)

    kernel = _build_corrected_kernel(nx_i, ny_i, dx_f, near_i, quad_i)
    fft_shape = (3 * nx_i - 2, 3 * ny_i - 2)
    kernel_fft = fft(_embed_kernel_for_linear_convolution(kernel, fft_shape))

    return CorrectedKernelDtN2D(nx_i, ny_i, dx_f, near_i, quad_i, kernel, kernel_fft, fft_shape)
end

"""
    build_periodic_dtn(nx, ny, dx)

Construct the exact periodic DtN operator using the Fourier symbol `|k|`.
"""
function build_periodic_dtn(nx::Integer, ny::Integer, dx::Real)
    nx_i = Int(nx)
    ny_i = Int(ny)
    dx_f = Float64(dx)

    symbol = _build_symbol(nx_i, ny_i, dx_f)
    kernel = _fftshift(real.(ifft(symbol)))

    return PeriodicDtN2D(nx_i, ny_i, dx_f, kernel, ComplexF64.(symbol))
end

"""Build the corrected real-space DtN kernel on a finite uniform grid."""
function _build_corrected_kernel(nx::Int, ny::Int, dx::Float64, near_radius::Int, quadrature_order::Int)
    kernel = zeros(Float64, 2 * nx - 1, 2 * ny - 1)
    cx = nx
    cy = ny
    alpha = _origin_laplacian_coefficient(dx)
    axis_weight = _axis_weight_with_quadratic_exactness(dx, quadrature_order, alpha)

    for n in -(ny - 1):(ny - 1), m in -(nx - 1):(nx - 1)
        if m == 0 && n == 0
            continue
        end

        weight = if (abs(m) == 1 && n == 0) || (m == 0 && abs(n) == 1)
            axis_weight
        else
            _offset_weight(m, n, dx, near_radius, quadrature_order)
        end
        kernel[cx + m, cy + n] = -weight
    end

    _add_origin_cell_correction!(kernel, cx, cy, alpha)
    kernel[cx, cy] = -sum(kernel) + kernel[cx, cy]

    return kernel
end

"""Return the discrete DtN weight for a given grid offset."""
function _offset_weight(m::Int, n::Int, dx::Float64, near_radius::Int, quadrature_order::Int)
    if max(abs(m), abs(n)) <= near_radius
        return _cell_average_weight(m, n, dx, quadrature_order)
    end

    r = dx * hypot(m, n)
    return DTN_CONSTANT * dx^2 / r^3
end

"""Compute a high-order cell-averaged singular-kernel weight for one offset cell."""
function _cell_average_weight(m::Int, n::Int, dx::Float64, quadrature_order::Int)
    nodes, weights = _gauss_legendre_rule(quadrature_order)
    acc = 0.0
    for (ξ, wξ) in zip(nodes, weights), (η, wη) in zip(nodes, weights)
        x = (m + 0.5 * ξ) * dx
        y = (n + 0.5 * η) * dx
        acc += wξ * wη / (x^2 + y^2)^(3 / 2)
    end
    return DTN_CONSTANT * (dx^2 / 4.0) * acc
end

"""Inject the analytic origin-cell Laplacian correction into the local kernel stencil."""
function _add_origin_cell_correction!(kernel::Matrix{Float64}, cx::Int, cy::Int, alpha::Float64)
    kernel[cx - 1, cy] += alpha
    kernel[cx + 1, cy] += alpha
    kernel[cx, cy - 1] += alpha
    kernel[cx, cy + 1] += alpha
    kernel[cx, cy] -= 4.0 * alpha
    return kernel
end

"""Return the origin-cell correction coefficient coming from the quadratic Taylor expansion."""
function _origin_laplacian_coefficient(dx::Float64)
    return -asinh(1.0) / (2.0 * π * dx)
end

"""Calibrate the nearest axis-neighbor weight to satisfy quadratic exactness locally."""
function _axis_weight_with_quadratic_exactness(dx::Float64, quadrature_order::Int, alpha::Float64)
    wdiag = _cell_average_weight(1, 1, dx, quadrature_order)
    i0 = 4.0 * dx * asinh(1.0)
    iaxis = _cell_average_inverse_radius(1, 0, dx, quadrature_order)
    idiag = _cell_average_inverse_radius(1, 1, dx, quadrature_order)
    rhs = -DTN_CONSTANT * (i0 + 4.0 * iaxis + 4.0 * idiag)
    return (-rhs + 4.0 * alpha * dx^2 - 8.0 * wdiag * dx^2) / (4.0 * dx^2)
end

"""Compute the cell average of `1/r` over one offset cell."""
function _cell_average_inverse_radius(m::Int, n::Int, dx::Float64, quadrature_order::Int)
    nodes, weights = _gauss_legendre_rule(quadrature_order)
    acc = 0.0
    for (ξ, wξ) in zip(nodes, weights), (η, wη) in zip(nodes, weights)
        x = (m + 0.5 * ξ) * dx
        y = (n + 0.5 * η) * dx
        acc += wξ * wη / sqrt(x^2 + y^2)
    end
    return (dx^2 / 4.0) * acc
end

"""Build the periodic DtN Fourier symbol `|k|` on the discrete grid."""
function _build_symbol(nx::Int, ny::Int, dx::Float64)
    lx = nx * dx
    ly = ny * dx
    kx = _angular_frequencies(nx, lx)
    ky = _angular_frequencies(ny, ly)

    symbol = zeros(Float64, nx, ny)
    for j in 1:ny, i in 1:nx
        symbol[i, j] = hypot(kx[i], ky[j])
    end
    return symbol
end

"""Return centered angular frequencies for a periodic grid of a given physical length."""
function _angular_frequencies(n::Int, length::Float64)
    freqs = zeros(Float64, n)
    half = fld(n, 2)
    scale = 2π / length
    for i in 1:n
        m = i - 1
        if m > half
            m -= n
        end
        freqs[i] = scale * m
    end
    return freqs
end

"""Apply a 2D FFT shift to move the zero mode to the center."""
function _fftshift(A::AbstractMatrix)
    sx = fld(size(A, 1), 2)
    sy = fld(size(A, 2), 2)
    return circshift(A, (sx, sy))
end

"""Embed a centered kernel into a padded array for linear FFT convolution."""
function _embed_kernel_for_linear_convolution(kernel::Matrix{Float64}, fft_shape::Tuple{Int, Int})
    padded = zeros(Float64, fft_shape...)
    nxk, nyk = size(kernel)
    padded[1:nxk, 1:nyk] .= kernel
    return padded
end

"""Return the supported Gauss-Legendre quadrature rule nodes and weights."""
function _gauss_legendre_rule(order::Int)
    if order == 8
        nodes = [
            -0.9602898564975363,
            -0.7966664774136267,
            -0.5255324099163290,
            -0.1834346424956498,
             0.1834346424956498,
             0.5255324099163290,
             0.7966664774136267,
             0.9602898564975363,
        ]
        weights = [
            0.1012285362903763,
            0.2223810344533745,
            0.3137066458778873,
            0.3626837833783620,
            0.3626837833783620,
            0.3137066458778873,
            0.2223810344533745,
            0.1012285362903763,
        ]
        return nodes, weights
    end

    throw(ArgumentError("unsupported quadrature order: $order"))
end

"""
    grid_coordinates(nx, ny, dx)

Return centered Cartesian grid coordinates as column/row arrays compatible with broadcasting.
"""
function grid_coordinates(nx::Integer, ny::Integer, dx::Real)
    xs = ((collect(0:Int(nx)-1) .- fld(Int(nx), 2)) .* Float64(dx))
    ys = ((collect(0:Int(ny)-1) .- fld(Int(ny), 2)) .* Float64(dx))
    return reshape(xs, Int(nx), 1), reshape(ys, 1, Int(ny))
end
