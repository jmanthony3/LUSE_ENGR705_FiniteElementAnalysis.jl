var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = FiniteElementAnalysis","category":"page"},{"location":"#FiniteElementAnalysis.jl","page":"Home","title":"FiniteElementAnalysis.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for FiniteElementAnalysis.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"#Functions","page":"Home","title":"Functions","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Modules = [FiniteElementAnalysis]","category":"page"},{"location":"#FiniteElementAnalysis.axialstresses-Tuple{Dict{String, FiniteElementAnalysis.LineElements}, Vector{Float64}}","page":"Home","title":"FiniteElementAnalysis.axialstresses","text":"axialstresses(elements::Dict{String, LineElements}, nodedisplacements::Vector{Float64}, dims=1)\naxialstresses(elements::Dict{String, BeamElements}, bendingmoments::Vector{Float64}, distancesfromneutralaxis::Union{<:Real, Vector{<:Real}}, dims=2)\n\nDetermine the axial stresses present in element-type system.\n\nIf line elements, of two-force members in their local coordinate systems. If beam elements, of beam elements which assumes no direct, axial forces.\n\nSee also: prepareelements_line, prepareelements_beam, solve, and transformationmatrixstar.\n\nArguments\n\nelements: the prepared dictionary according to element type.\nnodedisplacements: the vector for displacement of each node in the global coordinate system.\nbendingmoments: the vector of reaction moments at each node from solve.F.\ndistancesfromneutralaxis: the maximum distance(s) of nodes from neutral axis.\ndims=1: the number of dimensions being analyzed.\n\nExamples\n\nA, E, L = 1, 10e6, 100          # [in², psi, in]\nangles = [120, 0, 210]          # [°]\nL = L ./ abs.(cosd.(angles))    # [in]\nelements = Dict(\"1\"=>(1, 2), \"2\"=>(1, 3), \"3\"=>(1, 4))\npreparedelements = prepareelements_line(elements, A, E, L, angles, dims=2)\nnodeboundaryconditions = [2, 3, 4]\nF_applied = [(1, (1e3, 1e3))]   # (node, (force_x [lb], force_y [lb]))\nU = solve(preparedelements, nodeboundaryconditions, F_applied, dims=2).U\naxialstresses(preparedelements, U, dims=2)\n\nA, E, I, L = 0.004887968, 200e9, 102e6 * (1/1e3)^4, 6               # [m², Pa, m⁴, m]\nQ, c, t_w = 323473.474 * (1/1e3)^3, 176.5e-3, 6.48e-3               # [m³, m, m]\nW(x) = -10e3(x)                                                     # [N/m]\nF_applied = ([(W, 1, 2, 3, 4)])\nnumberofelements = 10\nL /= numberofelements\nx = 0:L:18\nelements, nodeboundaryconditions = Dict{String, Tuple{Int, Int}}(), []\nfor node ∈ 1:1:numberofelements+1\n    if node == 1\n        push!(nodeboundaryconditions, node)\n    else\n        elements[\"$(node-1)\"] = (node-1, node)\n        if x[node]%6 == 0\n            push!(nodeboundaryconditions, (node, (0, Inf)))\n        end\n    end\nend\npreparedelements = prepareelements_beam(elements, E, I, L)\nsol = solve(preparedelements, nodeboundaryconditions, F_applied)    # [m]\nV = [sol.F[2i - 1] for i ∈ 1:1:length(sol.F)÷2]                     # [N]\nM = [sol.F[2i] for i ∈ 1:1:length(sol.F)÷2]                         # [N-m]\nmaximum(axialstresses(preparedelements, M, c))                      # [Pa]\n\n\n\n\n\n","category":"method"},{"location":"#FiniteElementAnalysis.globalstiffnessmatrix-Tuple{Dict{String, FiniteElementAnalysis.LineElements}}","page":"Home","title":"FiniteElementAnalysis.globalstiffnessmatrix","text":"globalstiffnessmatrix(elements::Dict{String, LineElements}, dims::Integer=1)::Matrix{Float64}\nglobalstiffnessmatrix(elements::Dict{String, BeamElements}, dims::Integer=2)::Matrix{Float64}\n\nAssembles the global stiffness matrix from the dictionary, elements.\n\nSee also: prepareelements_line, prepareelements_beam, and transformationmatrix.\n\nArguments\n\nelements: entries infer localstiffnessmatrices and nodeconnections for each element in local coordinate system transformed by angle.\ndims: the number of dimensions being analyzed.\n\nExamples\n\nA = [0.75, 0.5, 1]   # [in²]\nE = 10.6e6           # [psi]\nL = [12, 9, 8]       # [in]\nelements = Dict(\"1\"=>(1, 2), \"2\"=>(2, 3), \"3\"=>(3, 4))\nglobalstiffnessmatrix(prepareelements_line(elements, A, E, L))\n\nE, I, L = 30e6, 510, 60 # [psi, in⁴, in]\nW(x) = -12\\1e3(x)      # [lb/in]\nF_applied = ([(W, 1, 2, 3)])\nelements = Dict(\"1\"=>(1, 2), \"2\"=>(2, 3))\nglobalstiffnessmatrix(prepareelements_beam(elements, E, I, L))\n\n\n\n\n\n","category":"method"},{"location":"#FiniteElementAnalysis.prepareelements_beam","page":"Home","title":"FiniteElementAnalysis.prepareelements_beam","text":"prepareelements_beam(elements, elementmoduli, elementareasofinertia, elementlengths, nodeangles=0, isdegrees=true; dims=2, stiffnesscoefficients=nothing)::Dict{String, BeamElements}\n\nReturn the dictionary with properties specific to those beam elements between node pairs.\n\nSee also: globalstiffnessmatrix, solve, and axialstresses.\n\nArguments\n\nelements::Dict{String, Tuple{Int, Int}}: the dictionary of \"element\"=>(node 1, node 2) for the element between node pairs.\nelementmoduli::Union{<:Real, Vector{<:Real}}: the elastic modulus of each element.\nelementareasofinertia::Union{<:Real, Vector{<:Real}}: the cross-sectional area of each element.\nelementlengths::Union{<:Real, Vector{<:Real}}: the length of each element.\nnodeangles::Union{<:Real, Vector{<:Real}}=0: the angle between the global and local coordinate systems of each element.\nisdegrees::Bool=true: control value whether elements of nodeangles is degrees (true) or radians (false).\ndims::Integer=2: the number of dimensions being analyzed.\nstiffnesscoefficients::Union{<:Real, Vector{<:Real}}: the stiffness coefficient, k associated with each element.\n\nExamples\n\nE, I, L = 30e6, 200, 20     # [psi, in⁴, ft]\nelements = Dict(\"1\"=>(1, 2),\"2\"=>(2, 3))\nprepareelements_beam(elements, E, I, L)\n\n\n\n\n\n","category":"function"},{"location":"#FiniteElementAnalysis.prepareelements_line","page":"Home","title":"FiniteElementAnalysis.prepareelements_line","text":"prepareelements_line(elements, elementareas, elementmoduli, elementlengths, nodeangles=0, isdegrees=true; dims=1, stiffnesscoefficients=nothing)::Dict{String, LineElements}\n\nReturn the dictionary with properties specific to those line elements between node pairs.\n\nSee also: globalstiffnessmatrix, solve, and axialstresses.\n\nArguments\n\nelements::Dict{String, Tuple{Int, Int}}: the dictionary of \"element\"=>(node 1, node 2) for the element between node pairs.\nelementareas::Union{<:Real, Vector{<:Real}}: the cross-sectional area(s) of elements.\nelementmoduli::Union{<:Real, Vector{<:Real}}: the elastic modul(us/i) of elements.\nelementlengths::Union{<:Real, Vector{<:Real}}: the length(s) of elements.\nnodeangles::Union{<:Real, Vector{<:Real}}=0: the angle between the global and local coordinate systems of elements.\nisdegrees::Bool=true: control value whether elements of nodeangles is degrees (true) or radians (false).\ndims::Integer=1: the number of dimensions being analyzed.\nstiffnesscoefficients::Union{<:Real, Vector{<:Real}}: the stiffness coefficient, k associated with each element.\n\nExamples\n\nA = [0.75, 0.5, 1]   # [in²]\nE = 10.6e6           # [psi]\nL = [12, 9, 8]       # [in]\nelements = Dict(\"1\"=>(1, 2), \"2\"=>(2, 3), \"3\"=>(3, 4))\nprepareelements_line(elements, A, E, L)\n\n\n\n\n\n","category":"function"},{"location":"#FiniteElementAnalysis.reducedglobalstiffnessmatrix-Tuple{Matrix{Float64}, Union{Vector{Integer}, Vector{Tuple{Integer, Integer}}}}","page":"Home","title":"FiniteElementAnalysis.reducedglobalstiffnessmatrix","text":"reducedglobalstiffnessmatrix(K::Matrix{Float64}, nodes::Union{Vector{Integer}, Vector{Tuple{Integer, Integer}}}, dims::Integer=1)::Matrix{Float64}\n\nEliminate the rows and columns of K according to nodes which implies the displacement of the appropriate node is known: e.g. is fixed or a prescribed displacement.\n\nSee also: globalstiffnessmatrix.\n\nExamples\n\nA = [0.75, 0.5, 1]   # [in²]\nE = 10.6e6           # [psi]\nL = [12, 9, 8]       # [in]\nelements = Dict(\"1\"=>(1, 2), \"2\"=>(2, 3), \"3\"=>(3, 4))\nK = globalstiffnessmatrix(prepareelements_line(elements, A, E, L))\nskipnodes = Vector{Integer}([1, 4])\nreducedglobalstiffnessmatrix(K, skipnodes)\n\n\n\n\n\n","category":"method"},{"location":"#FiniteElementAnalysis.shearstresses-Tuple{Dict{String, FiniteElementAnalysis.BeamElements}, Vector{Float64}, Union{Real, Vector{var\"#s80\"} where var\"#s80\"<:Real}, Union{Real, Vector{var\"#s78\"} where var\"#s78\"<:Real}}","page":"Home","title":"FiniteElementAnalysis.shearstresses","text":"shearstresses(elements, shearforces, firstmomentofareas, bendingthicknesses, dims=2)\n\nDetermine the axial stresses present in beam element system.\n\nSee also: prepareelements_beam, solve, and transformationmatrixstar.\n\nArguments\n\nelements::Dict{String, BeamElements}: the prepared dictionary according to element type.\nshearforces::Vector{Float64}: the vector for reaction, shear forces at each node in the global coordinate system.\nfirstmomentofarea::Union{<:Real, Vector{<:Real}}: the first moment of area(s), Q above the neutral axis.\nbendingthicknesses::Union{<:Real, Vector{<:Real}}: the web thickness of I-beam.\ndims::Integer=2: the number of dimensions being analyzed.\n\nExamples\n\nA, E, I, L = 0.004887968, 200e9, 102e6 * (1/1e3)^4, 6               # [m², Pa, m⁴, m]\nQ, c, t_w = 323473.474 * (1/1e3)^3, 176.5e-3, 6.48e-3               # [m³, m, m]\nW(x) = -10e3(x)                                                     # [N/m]\nF_applied = ([(W, 1, 2, 3, 4)])\nnumberofelements = 10\nL /= numberofelements\nx = 0:L:18\nelements, nodeboundaryconditions = Dict{String, Tuple{Int, Int}}(), []\nfor node ∈ 1:1:numberofelements+1\n    if node == 1\n        push!(nodeboundaryconditions, node)\n    else\n        elements[\"$(node-1)\"] = (node-1, node)\n        if x[node]%6 == 0\n            push!(nodeboundaryconditions, (node, (0, Inf)))\n        end\n    end\nend\npreparedelements = prepareelements_beam(elements, E, I, L)\nsol = solve(preparedelements, nodeboundaryconditions, F_applied)    # [m]\nV = [sol.F[2i - 1] for i ∈ 1:1:length(sol.F)÷2]                     # [N]\nM = [sol.F[2i] for i ∈ 1:1:length(sol.F)÷2]                         # [N-m]\nmaximum(shearstresses(preparedelements, V, Q, t_w))                 # [Pa]\n\n\n\n\n\n","category":"method"},{"location":"#FiniteElementAnalysis.solve-Tuple{Dict{String, FiniteElementAnalysis.LineElements}, Any, Union{Vector{var\"#s72\"} where var\"#s72\"<:(Tuple{Integer, var\"#s71\"} where var\"#s71\"<:Real), Vector{var\"#s70\"} where var\"#s70\"<:Tuple{Integer, Tuple{Real, Real}}}}","page":"Home","title":"FiniteElementAnalysis.solve","text":"solve(elements::Dict{String, LineElements}, nodeboundaryconditions, nodeforces::Vector{<:Tuple{Integer, Vararg{<:Real}}}, dims::Integer=1)\nsolve(elements::Dict{String, BeamElements}, nodeboundaryconditions, nodeforces::Vector{<:Tuple{Integer, Vararg{<:Real}}}, dims::Integer=2)\n\nPerform Finite Element Analysis (FEA) with entries of elements and applied boundary conditions by the Direct Stiffness Method according to Hooke's Law of linear-elastic deformation.\n\nIf line elements, stiffness coefficient defined by: k  AEL (typical of trusses with two-force members). If beam elements, stiffness coefficient defined by: k  EIL³.\n\nSee also: prepareelements_line, prepareelements_beam, globalstiffnessmatrix, and reducedglobalstiffnessmatrix.\n\nArguments\n\nelements: the prepared dictionary according to element type.\nnodeboundaryconditions: the tuple of nodes that either experience zero displacement (tuple item of integer for node) or prescribed displacement constrained by the stiffness of a certain element (tuple item of node integer and tuple of real displacement and integer element).\nnodeforces::Vector{<:Tuple{Integer, Vararg{<:Real}}}: the vector of tuple pairs for nodes experiencing some applied force.\ndims::Integer=1: the number of dimensions being analyzed.\nalg=nothing: algorithm of choice to solve linear system of equations according to LinearSolve.jl.\n\nExamples\n\n1D\n\nnodeboundaryconditions:= tuple of nodes with zero displacement.\n\nA, E, L = 4, 30e6, 30               # [in², psi, in]\nelements = Dict(\"1\"=>(1, 2), \"2\"=>(2, 3), \"3\"=>(3, 4))\npreparedelements = prepareelements_line(elements, A, E, L)\nnodeboundaryconditions = (1, 4)\nF_applied = [(2, 5e3), (3, -10e3)]  # (node, force [lb])\nsolve(preparedelements, nodeboundaryconditions, F_applied).U\n\nnodeboundaryconditions:= tuple of nodes with zero and prescribed displacements.\n\nThe tuple items are identified in the following way:\n\nNode (1) with zero displacement.\nNode (3) has a prescribed displacement of 25e-3 constrained by the deformation of element 2.\n\nA, E, L = 4e-4, 210e9, 2                    # [m², Pa, m]\nelements = Dict(\"1\"=>(1, 2), \"2\"=>(2, 3))\npreparedelements = prepareelements_line(elements, A, E, L)\nnodeboundaryconditions = [1, (3, 25e-3, 2)] # (node, displacement [m], element)\nF_applied = [(2, -5e3)]                     # (node, force [N])\nsolve(preparedelements, nodeboundaryconditions, F_applied).U\n\n2D\n\nA, E, L = 1, 10e6, 100          # [in², psi, in]\nangles = [120, 0, 210]          # [°]\nL = L ./ abs.(cosd.(angles))    # [in]\nelements = Dict(\"1\"=>(1, 2), \"2\"=>(1, 3), \"3\"=>(1, 4))\npreparedelements = prepareelements_line(elements, A, E, L, angles, dims=2)\nnodeboundaryconditions = [2, 3, 4]\nF_applied = [(1, (1e3, 1e3))]   # (node, (force_x [lb], force_y [lb]))\nsolve(preparedelements, nodeboundaryconditions, F_applied, dims=2).U\n\n\n\n\n\n","category":"method"},{"location":"#Index","page":"Home","title":"Index","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"}]
}
