using Test
using VerticalDynamics
using LinearAlgebra

include("test_api.jl")
include("test_kernel_structure.jl")
include("test_impulse_response.jl")
include("test_fft_apply.jl")
include("test_convergence.jl")
include("test_visualization.jl")
include("test_solver_api.jl")
include("test_solver_domain.jl")
include("test_solver_assembly.jl")
include("test_solver_operator.jl")
include("test_solver_step.jl")
