@testset "Visualization helpers" begin
    history = [
        SolverState(0.0, [0.0 1.0; -2.0 0.5], zeros(2, 2), zeros(2, 2), zeros(3), zeros(3)),
        SolverState(0.1, [3.0 -1.0; 0.2 0.0], zeros(2, 2), zeros(2, 2), zeros(3), zeros(3)),
    ]

    @test compute_eta_color_limits(history) == (-2.0, 3.0)
end
