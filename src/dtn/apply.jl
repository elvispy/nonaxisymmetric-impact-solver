"""
    apply(op, x)

Apply a DtN operator. The corrected nonperiodic operator uses zero-padded linear convolution.
"""
function apply_direct(op::CorrectedKernelDtN2D, x::AbstractMatrix{<:Real})
    size(x) == (op.nx, op.ny) || throw(DimensionMismatch("input shape must match operator grid"))
    y = zeros(Float64, op.nx, op.ny)
    cx = op.nx
    cy = op.ny

    for j in 1:op.ny, i in 1:op.nx
        acc = 0.0
        for q in 1:op.ny, p in 1:op.nx
            acc += op.kernel[cx + (i - p), cy + (j - q)] * Float64(x[p, q])
        end
        y[i, j] = acc
    end

    return y
end

"""Apply the nonperiodic corrected DtN operator by FFT-based linear convolution."""
function apply(op::CorrectedKernelDtN2D, x::AbstractMatrix{<:Real})
    size(x) == (op.nx, op.ny) || throw(DimensionMismatch("input shape must match operator grid"))
    padded_x = zeros(Float64, op.fft_shape...)
    padded_x[1:op.nx, 1:op.ny] .= Float64.(x)

    y_full = real.(ifft(op.kernel_fft .* fft(padded_x)))
    i0 = op.nx
    j0 = op.ny
    return y_full[i0:(i0 + op.nx - 1), j0:(j0 + op.ny - 1)]
end
