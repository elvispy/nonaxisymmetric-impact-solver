using Plots
using Printf

"""
    compute_eta_color_limits(history)

Return a consistent `(zmin, zmax)` pair for heatmap rendering across a solver history.
"""
function compute_eta_color_limits(history::AbstractVector{<:SolverState})
    isempty(history) && throw(ArgumentError("history must contain at least one state"))
    zmin = minimum(minimum(state.eta) for state in history)
    zmax = maximum(maximum(state.eta) for state in history)

    if zmin == zmax
        pad = zmin == 0.0 ? 1.0 : 0.05 * abs(zmin)
        return zmin - pad, zmax + pad
    end

    return zmin, zmax
end

"""
    render_eta_video(history, x, y, output_path; fps=20, title_prefix="eta", frames_dir=nothing)

Render top-view heatmap frames of `η` and assemble them into an MP4 with `ffmpeg`.
"""
function render_eta_video(
    history::AbstractVector{<:SolverState},
    x::AbstractMatrix{<:Real},
    y::AbstractMatrix{<:Real},
    output_path::AbstractString;
    fps::Integer=20,
    title_prefix::AbstractString="eta",
    frames_dir::Union{Nothing, AbstractString}=nothing,
)
    isempty(history) && throw(ArgumentError("history must contain at least one state"))
    size(first(history).eta) == size(x) == size(y) || throw(DimensionMismatch("history fields and grid must match"))

    ffmpeg = Sys.which("ffmpeg")
    isnothing(ffmpeg) && throw(ArgumentError("ffmpeg was not found in PATH"))

    output_dir = dirname(output_path)
    isempty(output_dir) || mkpath(output_dir)

    frame_root = isnothing(frames_dir) ? joinpath(output_dir, "$(splitext(basename(output_path))[1])_frames") : String(frames_dir)
    mkpath(frame_root)

    xs = vec(x[:, 1])
    ys = vec(y[1, :])
    clim = compute_eta_color_limits(history)

    default(size=(800, 650))
    for (idx, state) in enumerate(history)
        frame_path = joinpath(frame_root, @sprintf("frame_%05d.png", idx - 1))
        plt = heatmap(
            xs,
            ys,
            transpose(state.eta);
            clim=clim,
            color=:balance,
            xlabel="x",
            ylabel="y",
            aspect_ratio=:equal,
            title="$(title_prefix), t=$(round(state.time, digits=4))",
            colorbar_title="η",
        )
        savefig(plt, frame_path)
    end

    cmd = `$(ffmpeg) -y -framerate $(Int(fps)) -i $(joinpath(frame_root, "frame_%05d.png")) -pix_fmt yuv420p -vcodec libx264 $(output_path)`
    run(cmd)
    return output_path
end
