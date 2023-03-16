using FiniteElementAnalysis
using Test

@testset "FiniteElementAnalysis.jl" begin
    # ENGR 705-001: Finite Element Analysis - Homework 3
    ## 3.1
    A, E, L = 4, 30e6, 30 # [in², psi, in]
    elements = Dict("1"=>(1, 2), "2"=>(2, 3), "3"=>(3, 4))
    nodeboundaryconditions = (1, 4)
    F_applied = [(2, 5e3), (3, -10e3)]
    preparedelements = prepareelements_line(elements, A, E, L)
    U = solve(preparedelements, nodeboundaryconditions, F_applied).U
    @test round.(U, digits=6) == [0, 0, -0.00125, 0]
    ## 3.2
    A, E, L = 4e-4, 210e9, 2 # [m², Pa, m]
    elements = Dict("1"=>(1, 2), "2"=>(2, 3))
    nodeboundaryconditions = [1, (3, 25e-3, 2)]
    F_applied = [(2, -5e3)]
    preparedelements = prepareelements_line(elements, A, E, L)
    U = solve(preparedelements, nodeboundaryconditions, F_applied).U
    @test round.(U, digits=6) == [0, 0.012440, 0.025]
    ## 3.3
    A, E, L = 1, 10e6, 100          # [in², psi, in]
    angles = [120, 0, 210]          # [°]
    L = L ./ abs.(cosd.(angles))
    elements = Dict("1"=>(1, 2), "2"=>(1, 3), "3"=>(1, 4))
    nodeboundaryconditions = [2, 3, 4]
    F_applied = [(1, (1e3, 1e3))]   # (node, (force_x [lb], force_y [lb]))
    preparedelements = prepareelements_line(elements, A, E, L, angles, dims=2)
    sol = solve(preparedelements, nodeboundaryconditions, F_applied, dims=2)
    σ = axialstresses(preparedelements, sol.U, dims=2)
    @test round.(σ, digits=6) == [-577.350269, -422.649731, 1000]

    # ENGR 705-001: Finite Element Analysis - FEA 1
    ## FEA1.1
    A = [0.75, 0.5, 1]                  # [in²]
    E = 10.6e6                          # [psi]
    L = [12, 9, 8]                      # [in]
    elements = Dict("1"=>(1, 2), "2"=>(2, 3), "3"=>(3, 4))
    nodeboundaryconditions = (1, 4)
    F_applied = [(2, 4e3), (3, 1.5e3)]  # [lb]
    preparedelements = prepareelements_line(elements, A, E, L)
    U = solve(preparedelements, nodeboundaryconditions, F_applied).U
    # @test round(U[2], digits=6) == 4.169e-03
    σ = axialstresses(preparedelements, U)
    @test round.(σ, digits=3) == [3682.540, -2476.190, -2738.095]
    ## FEA1.2
    A = [0.75, 2.15, 2.95]                          # [in²]
    E = 29e6                                        # [psi]
    L = [5, 5, 10] .* 12                            # [in]
    elements = Dict("1"=>(1, 2), "2"=>(2, 3), "3"=>(3, 4))
    nodeboundaryconditions = (4)
    F_applied = [(1, -100), (2, -150), (3, -200)]   # [lb]
    preparedelements = prepareelements_line(elements, A, E, L)
    U = solve(preparedelements, nodeboundaryconditions, F_applied).U
    # @test round.(U, digits=6) == [-1.148e-03, -8.718e-04, -6.312e-04, 0]
    σ = axialstresses(preparedelements, U)
    @test round.(σ, digits=3) == [133.333, 116.279, 152.542]
    ## FEA1.3
    A, E, L, θ = 2, 29e6, 8*12, 60          # [in², psi, in, °]
    angles = [0, 0, 180-θ, θ, 180-θ, θ, 0]  # [°]
    L = L ./ cosd.(angles)                  # [in]
    elements = Dict(
        "1"=>(1, 2),
        "2"=>(2, 3),
        "3"=>(3, 5),
        "4"=>(2, 5),
        "5"=>(2, 4),
        "6"=>(1, 4),
        "7"=>(4, 5)
    )
    nodeboundaryconditions = [1, (3, (Inf, 0))]
    F_applied = [(2, (0, -1e3))]            # (node, (force_x [lb], force_y [lb]))
    preparedelements = prepareelements_line(elements, A, E, L, angles, dims=2)
    U = solve(preparedelements, nodeboundaryconditions, F_applied, dims=2).U
    # U_x = [U[2i - 1] for i ∈ 1:1:length(U)÷2]
    # U_y = [U[2i] for i ∈ 1:1:length(U)÷2]
    # @test round.(U_x, digits=6) == [0, 4.778e-04, 9.556e-04, 9.556e-04, 0]
    # @test round.(U_y, digits=6) == [0, -8.276e-04, 0, -1.103e-03, 0]
    σ = axialstresses(preparedelements, U, dims=2)
    @test round.(σ, digits=3) == [144.338, 144.338, -288.675, 288.675, 288.675, -288.675, -288.675]
    ## FEA1.4
    A, E = 10, 29e6                 # [in², psi]
    L = [                           # [ft]
        15,                         # (1 -> 2)
        15,                         # (1 -> 3)
        15√(2),                     # (1 -> 4)
        15,                         # (2 -> 4)
        15,                         # (3 -> 4)
        √(2.5^2 + 10^2),            # (3 -> 5)
        √(12.5^2 + 10^2),           # (4 -> 5)
        √(2.5^2 + 10^2),            # (4 -> 6)
        10,                         # (5 -> 6)
        10,                         # (5 -> 7)
        10√(2),                     # (5 -> 8)
        10,                         # (6 -> 8)
        10,                         # (7 -> 8)
        √(5^2 + 10^2),              # (7 -> 9)
        √(5^2 + 10^2),              # (9 -> 11)
        5,                          # (9 -> 12)
        √(5^2 + 10^2),              # (9 -> 13)
        10,                         # (11 -> 12)
        10,                         # (12 -> 13)
        10,                         # (7 -> 13)
        √(5^2 + 10^2),              # (8 -> 10)
        10√(2),                     # (8 -> 13)
        10,                         # (8 -> 14)
        10,                         # (13 -> 14)
        10,                         # (10 -> 14)
        5,                          # (10 -> 15)
        10,                         # (14 -> 15)
        √(5^2 + 10^2),              # (10 -> 16)
        10,                         # (15 -> 16)
    ]
    theta(v1, v2) = acosd((v1*v2')/(√(sum(v1 .^ 2))*√(sum(v2 .^ 2))))[1]
    angles = [                      # [°]
        0,                          # (1 -> 2)
        90,                         # (1 -> 3)
        45,                         # (1 -> 4)
        90,                         # (2 -> 4)
        0,                          # (3 -> 4)
        theta([1 0], [2.5 10]),     # (3 -> 5)
        theta([1 0], [-12.5 10]),   # (4 -> 5)
        theta([1 0], [-2.5 10]),    # (4 -> 6)
        0,                          # (5 -> 6)
        90,                         # (5 -> 7)
        45,                         # (5 -> 8)
        90,                         # (6 -> 8)
        0,                          # (7 -> 8)
        theta([1 0], [-10 5]),      # (7 -> 9)
        theta([1 0], [-10 5]),      # (9 -> 11)
        90,                         # (9 -> 12)
        theta([1 0], [10 5]),       # (9 -> 13)
        0,                          # (11 -> 12)
        0,                          # (12 -> 13)
        90,                         # (7 -> 13)
        theta([1 0], [10 5]),       # (8 -> 10)
        theta([1 0], [-10 10]),     # (8 -> 13)
        90,                         # (8 -> 14)
        0,                          # (13 -> 14)
        theta([1 0], [-10 5]),      # (10 -> 14)
        90,                         # (10 -> 15)
        0,                          # (14 -> 15)
        theta([1 0], [10 5]),       # (10 -> 16)
        0,                          # (15 -> 16)
    ]
    elements = Dict(
        "1"=>(1, 2),
        "2"=>(1, 3),
        "3"=>(1, 4),
        "4"=>(2, 4),
        "5"=>(3, 4),
        "6"=>(3, 5),
        "7"=>(4, 5),
        "8"=>(4, 6),
        "9"=>(5, 6),
        "10"=>(5, 7),
        "11"=>(5, 8),
        "12"=>(6, 8),
        "13"=>(7, 8),
        "14"=>(7, 9),
        "15"=>(9, 11),
        "16"=>(9, 12),
        "17"=>(9, 13),
        "18"=>(11, 12),
        "19"=>(12, 13),
        "20"=>(7, 13),
        "21"=>(8, 10),
        "22"=>(8, 13),
        "23"=>(8, 14),
        "24"=>(13, 14),
        "25"=>(10, 14),
        "26"=>(10, 15),
        "27"=>(14, 15),
        "28"=>(10, 16),
        "29"=>(15, 16),
    )
    nodeboundaryconditions = [(1, (Inf, 0)), 2]
    F_applied = [                   # (node, (force_x [lb], force_y [lb]))
        (11, (0, -1e3)),
        (16, (0, -1e3))
    ]
    preparedelements = prepareelements_line(elements, A, E, L, angles, dims=2)
    U = solve(preparedelements, nodeboundaryconditions, F_applied, dims=2).U
    # U_x = [U[2i - 1] for i ∈ 1:1:length(U)÷2]
    # U_y = [U[2i] for i ∈ 1:1:length(U)÷2]
    # @test round.(U_x, digits=6) == [0, 0, 3.233e-06, 4.310e-06, 2.033e-06, 1.314e-06, 1.170e-05, 5.955e-06, 1.337e-05, 5.337e-06, -4.490e-06, 1.257e-06, 7.004e-06, 1.275e-05, 1.850e-05, 2.425e-05]
    # @test round.(U_y, digits=6) == [0, 0, -4.310e-06, -4.310e-06, -7.157e-06, -8.206e-06, -1.003e-05, -1.108e-05, -2.276e-05, -2.591e-05, -7.454e-05, -2.276e-05, -1.003e-05, -1.108e-05, -2.591e-05, -7.979e-05]
    σ = axialstresses(preparedelements, U, dims=2)
    @test round.(σ, digits=3) == [0, -100.000, 0, -100.000, -200.000, -223.607, -223.607, 0, 0, 200.000, 200.000, -100.000, 0, -223.607, 0, 0, 200.000, 0, 0, 200.000, -223.607, 200.000, 0, -100.000, 25.000, -103.078, 0, -103.078, -25.000]

    # ENGR 705-001: Finite Element Analysis - Set 5-3
    ## Set 5-3.1
    E, I = 30e6, 510            # [psi, in⁴]
    numberofelements = 2
    L = 120/numberofelements    # [in]
    x = 0:5:10
    W(x) = -12\1.e3(x)          # [lb/in]
    F_applied = [(W, 1, 2, 3)]  # (node, (force_x [lb], force_y [lb]))
    elements = Dict{String, Tuple{Int, Int}}()
    nodeboundaryconditions = []
    for node ∈ 1:1:numberofelements+1
        if node == 1
            push!(nodeboundaryconditions, node)
        else
            elements["$(node-1)"] = (node-1, node)
            if x[node]%L == 0
                push!(nodeboundaryconditions, (node, (Inf, 0)))
            end
        end
    end
    preparedelements = prepareelements_beam(elements, E, I, L)
    U = solve(preparedelements, nodeboundaryconditions, F_applied).U
    @test round.(U, digits=6) == [0, 0, -0.05, -0.001373, -0.141176, -0.001569]
end
